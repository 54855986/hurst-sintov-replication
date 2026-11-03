# =============================================================================
# Power analysis by simulation
# Being made to feel guilty, anticipated pleasantness, moral rightness and donation
# =============================================================================
#
# PURPOSE
# This script estimates, for every hypothesis (H1 to H9), the probability of
#   (a) detecting an effect of the smallest size of interest, d = 0.20, and
#   (b) showing that an effect is absent when its true size is 0,
# for four total sample sizes and four ways of allocating participants
# between the two parts of the design. It then reports, for each sample size,
# the allocation that gives the lowest-powered hypothesis the highest power.
#
# DESIGN
# Four groups, assigned at random:
#   guilt induction (made to feel guilty vs. control)
#   x anticipation measures (not answered vs. answered before the choice).
# Groups that do not answer the anticipation measures complete the original
# survey of Hurst and Sintov (2022, Study 2); they carry H2 to H4.
# Groups that answer them rate, for keeping and for donating a possible $25
# bonus, how pleasant each option would feel and how right it would be;
# they carry H5 to H9. H1 uses all participants.
#
# HYPOTHESES AND MODELS
#   H1  induction -> felt guilt (+)                    linear regression, all participants
#   H2  induction -> donation (-)                      logistic regression, original survey
#   H3  felt guilt -> donation (+)                     logistic regression, original survey
#   H4  induction -> felt guilt -> donation (+)        indirect effect, original survey
#   H5  induction -> pleasantness difference (+)       linear regression, anticipation measures
#   H6  induction -> rightness difference (+)          linear regression, anticipation measures
#   H7  pleasantness difference -> donation (-)        logistic regression, anticipation measures
#   H8  rightness difference -> donation (+)           same model
#   H9  H7 is stronger than H8: b_H7 + b_H8 < 0        same model
#   Pleasantness difference = keeping minus donating; rightness difference =
#   donating minus keeping; both standardized.
#
# DECISION RULES (95% confidence intervals)
#   Detected  the 95% CI does not include 0 and the effect points in the predicted direction.
#   Absent    the 95% CI does not include an effect of d = 0.20 in the predicted direction.
#   H4: detected if the 95% Monte Carlo CI of the indirect effect lies above 0;
#       absent if the 95% CI of either of its two links does not include d = 0.20.
#
# ASSUMPTIONS AND SIMPLIFICATIONS
#   * Control donation rate 43%, the rate in the control group of Hurst and Sintov (2022, Study 2).
#   * Environmental self-identity is the only covariate. The planned analysis also includes
#     gender, age, political orientation and education, which change power very little.
#   * The correlation between the pleasantness and the rightness difference is set by r_AM.
#     No estimate exists yet; 0 is the default. Rerun with another value to check H7 to H9.
#   * The indirect effect (H4) uses a Monte Carlo CI instead of the bootstrap of the analysis.
#   * Log odds are converted to d with d = log OR x sqrt(3) / pi (Chinn, 2000).
#
# HOW TO RUN
#   Set n_sim (1,000 for a first look, 5,000 for the reported results) and run the whole
#   script. Results are printed and saved next to the script as power_results.csv and
#   power_allocation_summary.csv. Running time grows with n_sim; with 5,000 it can take
#   an hour or more.
# =============================================================================

# ---- Settings -------------------------------------------------------------------
n_sim       <- 5000                         # simulated data sets per scenario
N_total     <- c(2000, 2200, 2400, 2600)    # total sample sizes
share_orig  <- c(600, 650, 700, 750) / 1200 # share of each condition in the original-survey groups
                                            # (at N = 2,400: 600/600, 650/550, 700/500, 750/450 per group)
sesoi       <- 0.20                         # smallest effect size of interest (d)
p_control   <- 0.43                         # donation rate in the control group
r_AM        <- 0                            # correlation of the pleasantness and rightness differences
n_mc        <- 2000                         # Monte Carlo draws for the H4 interval
set.seed(2026)

# ---- Helpers ---------------------------------------------------------------------
z95  <- qnorm(0.975)
to_logodds <- function(d) d * pi / sqrt(3)
to_d       <- function(b) b * sqrt(3) / pi

# Decision for one estimate (d units): direction = +1 or -1
decide <- function(est, se, direction) {
  lo <- unname(est - z95 * se); hi <- unname(est + z95 * se)
  if (direction > 0) c(detected = lo > 0, absent = hi < sesoi)
  else               c(detected = hi < 0, absent = lo > -sesoi)
}

# ---- One simulated data set --------------------------------------------------------
# eff: true effects in d units
simulate_data <- function(n_orig, n_ant, eff) {
  induction <- c(rep(0:1, each = n_orig), rep(0:1, each = n_ant))
  anticip   <- rep(0:1, times = c(2 * n_orig, 2 * n_ant))
  n   <- length(induction)
  eid <- rnorm(n)
  guilt <- eff$h1 * induction + 0.3 * eid + rnorm(n)
  pride <- -0.5 * induction + 0.3 * eid + rnorm(n)
  e1 <- rnorm(n); e2 <- rnorm(n)
  dA <- eff$h5 * induction + e1                                    # pleasantness difference
  dM <- eff$h6 * induction + 0.2 * eid + r_AM * e1 + sqrt(1 - r_AM^2) * e2  # rightness difference
  b_guilt  <- to_logodds(eff$h3)
  c_direct <- to_logodds(eff$h2) - eff$h1 * b_guilt                # total effect of the induction = h2
  logit <- ifelse(anticip == 0,
                  qlogis(p_control) + c_direct * induction + b_guilt * guilt + 0.3 * eid,
                  qlogis(p_control) + to_logodds(eff$h7) * dA + to_logodds(eff$h8) * dM + 0.3 * eid)
  data.frame(induction, anticip, eid, guilt, pride, dA, dM, donate = rbinom(n, 1, plogis(logit)))
}

# ---- All tests on one data set -----------------------------------------------------
run_tests <- function(dat) {
  z <- function(x) as.numeric(scale(x))
  orig <- dat[dat$anticip == 0, ]; ant <- dat[dat$anticip == 1, ]
  orig$guilt <- z(orig$guilt); orig$pride <- z(orig$pride)
  ant$dA <- z(ant$dA); ant$dM <- z(ant$dM)
  se <- function(m, term) sqrt(vcov(m)[term, term])

  m1 <- lm(z(guilt) ~ induction + eid, data = dat)
  m2 <- glm(donate ~ induction + eid, family = binomial, data = orig)
  ma <- lm(guilt ~ induction + eid, data = orig)
  mb <- glm(donate ~ induction + guilt + pride + eid, family = binomial, data = orig)
  m5 <- lm(dA ~ induction + eid, data = ant)
  m6 <- lm(dM ~ induction + eid, data = ant)
  m7 <- glm(donate ~ induction + dA + dM + eid, family = binomial, data = ant)

  a <- coef(ma)[["induction"]]; se_a <- se(ma, "induction")
  b <- to_d(coef(mb)[["guilt"]]); se_b <- to_d(se(mb, "guilt"))
  ab <- quantile(rnorm(n_mc, a, se_a) * rnorm(n_mc, b, se_b), c(0.025, 0.975))
  V <- vcov(m7)

  rbind(
    H1 = decide(coef(m1)[["induction"]], se(m1, "induction"), +1),
    H2 = decide(to_d(coef(m2)[["induction"]]), to_d(se(m2, "induction")), -1),
    H3 = decide(b, se_b, +1),
    H4 = c(detected = unname(ab[1] > 0),
           absent   = (a + z95 * se_a < sesoi) || (b + z95 * se_b < sesoi)),
    H5 = decide(coef(m5)[["induction"]], se(m5, "induction"), +1),
    H6 = decide(coef(m6)[["induction"]], se(m6, "induction"), +1),
    H7 = decide(to_d(coef(m7)[["dA"]]), to_d(se(m7, "dA")), -1),
    H8 = decide(to_d(coef(m7)[["dM"]]), to_d(se(m7, "dM")), +1),
    H9 = decide(to_d(coef(m7)[["dA"]] + coef(m7)[["dM"]]),
                to_d(sqrt(V["dA", "dA"] + V["dM", "dM"] + 2 * V["dA", "dM"])), -1))
}

# ---- Scenarios -----------------------------------------------------------------------
# effects: every tested effect equals d = 0.20 -> power to DETECT H1-H8;
#          H7 = -0.20 and H8 = +0.20 are equally strong -> power to show H9 ABSENT
# h9:      H7 = -0.40, H8 = +0.20 (H7 stronger by 0.20) -> power to DETECT H9
# null:    every tested effect is 0 -> power to show H1-H8 ABSENT
scenarios <- list(
  effects = list(h1 = 0.2, h2 = -0.2, h3 = 0.2, h5 = 0.2, h6 = 0.2, h7 = -0.2, h8 = 0.2),
  h9      = list(h1 = 0.2, h2 = -0.2, h3 = 0.2, h5 = 0.2, h6 = 0.2, h7 = -0.4, h8 = 0.2),
  null    = list(h1 = 0,   h2 = 0,    h3 = 0,   h5 = 0,   h6 = 0,   h7 = 0,    h8 = 0))

power_for <- function(n_orig, n_ant) {
  share <- lapply(scenarios, function(eff)
    apply(replicate(n_sim, run_tests(simulate_data(n_orig, n_ant, eff))), c(1, 2), mean))
  H <- rownames(share$effects)
  data.frame(H,
             detect = ifelse(H == "H9", share$h9[, "detected"], share$effects[, "detected"]),
             absent = ifelse(H == "H9", share$effects[, "absent"], share$null[, "absent"]))
}

# ---- Run: every sample size x every allocation -----------------------------------------
results <- list()
for (N in N_total) for (s in share_orig) {
  n_orig <- round(N / 2 * s); n_ant <- N / 2 - n_orig
  cat(sprintf("N = %d: %d per original-survey group, %d per anticipation group ...\n", N, n_orig, n_ant))
  results[[paste(N, n_orig)]] <- cbind(N = N, n_orig = n_orig, n_ant = n_ant, power_for(n_orig, n_ant))
}
res <- do.call(rbind, results)
res$detect <- round(res$detect, 3); res$absent <- round(res$absent, 3)

# ---- Summary ------------------------------------------------------------------------------
# Lowest power of each design (across H1-H9 and both kinds of power)
low <- aggregate(cbind(detect, absent) ~ N + n_orig + n_ant, data = res, FUN = min)
low$lowest <- pmin(low$detect, low$absent)
low <- low[order(low$N, low$n_orig), ]
# Best allocation per sample size: the one with the highest lowest power
best <- do.call(rbind, lapply(split(low, low$N), function(x) x[which.max(x$lowest), ]))

cat("\nPower at d = 0.20 with 95% CIs (", n_sim, " simulated data sets per scenario)\n\n", sep = "")
print(res, row.names = FALSE)
cat("\nLowest power of each design:\n"); print(low, row.names = FALSE)
cat("\nBest allocation for each sample size:\n"); print(best, row.names = FALSE)

# ---- Save -----------------------------------------------------------------------------------
out_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  dirname(rstudioapi::getSourceEditorContext()$path) else getwd()
write.csv(res, file.path(out_dir, "power_results.csv"), row.names = FALSE)
write.csv(low, file.path(out_dir, "power_allocation_summary.csv"), row.names = FALSE)
cat("\nSaved power_results.csv and power_allocation_summary.csv in", out_dir, "\n")
sessionInfo()
