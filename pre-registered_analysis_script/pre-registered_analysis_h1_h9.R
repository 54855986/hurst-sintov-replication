# =============================================================================
# Preregistered analysis of H1 to H9
# Being made to feel guilty, anticipated pleasantness, moral rightness and donation
# =============================================================================
#
# PURPOSE
# Runs the five models that test the nine hypotheses and gives each hypothesis
# one of four outcomes: supported, supported but too small to matter,
# not supported, or inconclusive.
#
# DESIGN
# Four groups, assigned at random:
#   induction  1 = made to feel guilty (negative feedback of Hurst and Sintov, 2022), 0 = control
#   anticip    1 = rated keeping and donating before the choice, 0 = completed the original survey only
# Groups with anticip = 0 carry H2 to H4; groups with anticip = 1 carry H5 to H9;
# H1 uses all participants.
#
# FIVE MODELS
#   Model 1  H1      all participants   linear:   felt guilt ~ induction
#   Model 2  H2      anticip = 0        logistic: donation ~ induction
#   Model 3  H3, H4  anticip = 0        PROCESS model 4: induction -> felt guilt and pride (parallel) -> donation
#   Model 4  H5, H6  anticip = 1        linear:   pleasantness difference ~ induction; rightness difference ~ induction
#   Model 5  H7-H9   anticip = 1        logistic: donation ~ induction + pleasantness difference + rightness difference
# All models include environmental self-identity, gender, age, political orientation
# and education as covariates, as in Hurst and Sintov (2022).
#
# DECISION RULES (95% confidence intervals, smallest effect size of interest d = 0.20)
#   Supported      the 95% CI does not include 0 and the effect points in the predicted direction.
#   Not supported  the 95% CI does not include an effect of d = 0.20 in the predicted direction.
#   Supported but too small to matter   both of the above.
#   Inconclusive   neither.
#   H4 uses the 95% bootstrap CI of the indirect effect; it is not supported if the
#   95% CI of either of its two links does not include d = 0.20.
#   H9 tests whether the sum of the H7 and H8 effects is below 0, that is, whether the
#   pleasantness difference weighs more than the rightness difference. It can be read
#   only if H7 and H8 point in their predicted directions.
#
# EFFECT SIZES
#   Continuous variables are standardized, so linear effects are Cohen's d. Logistic
#   effects are converted with effectsize::oddsratio_to_d(log = TRUE).
#
# HOW TO RUN
#   1. Save PROCESS for R as 05_analysis/PROCESS/process.R (download from processmacro.org).
#   2. Install effectsize once: install.packages("effectsize").
#   3. Part 1 creates one simulated data set so the code can be checked before data
#      collection. For the real data, replace Part 1 by reading the survey export.
#   4. Run the whole script. Results are printed and saved as analysis_results.csv.
# =============================================================================

library(effectsize)

n_boot    <- 5000                                                   # PROCESS bootstrap samples
boot_seed <- 2000000                                                # PROCESS seed
sesoi     <- 0.20                                                   # smallest effect size of interest (d)
covs      <- c("eid", "gender", "age", "politics", "education")     # covariates
set.seed(2026)

find_root <- function() {
  starts <- getwd()
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    starts <- c(dirname(rstudioapi::getSourceEditorContext()$path), starts)
  for (dir in starts) {
    while (!dir.exists(file.path(dir, "05_analysis")) && dirname(dir) != dir) dir <- dirname(dir)
    if (dir.exists(file.path(dir, "05_analysis"))) return(dir)
  }
  stop("Project folder not found. Set the working directory to the project folder.")
}
root <- find_root()
source(file.path(root, "05_analysis/PROCESS/process.R"))

# =============================================================================
# PART 1. DATA
# =============================================================================
# Simulated data with the planned structure. The effects are set by hand and are
# NOT predictions. For the real data, replace this part with:
#   dat <- read.csv("path/to/survey_export.csv")
# keep participants who passed the bot check and have a unique Prolific ID, and
# make sure these variables exist:
#   induction, anticip   0/1
#   guilt, pride         mean of three items each (1-5)
#   A_keep, A_donate     how pleasant keeping / donating would feel (0-100), anticip = 1 only
#   M_keep, M_donate     how right keeping / donating would be (1-7), anticip = 1 only
#   donate               1 = donate, 0 = keep
#   eid                  environmental self-identity, mean of three items (1-7)
#   gender               1 = woman, 0 = otherwise
#   age                  years
#   politics             political orientation (1 = very liberal to 7 = very conservative)
#   education            highest education (1-6)

n_orig <- 700; n_ant <- 500                                         # participants per group
dat <- data.frame(induction = c(rep(0:1, each = n_orig), rep(0:1, each = n_ant)),
                  anticip   = rep(0:1, times = c(2 * n_orig, 2 * n_ant)))
n <- nrow(dat)
clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)
dat$eid       <- clip(round(rnorm(n, 5, 1.2), 1), 1, 7)
dat$gender    <- rbinom(n, 1, 0.5)
dat$age       <- round(runif(n, 18, 75))
dat$politics  <- sample(1:7, n, replace = TRUE)
dat$education <- sample(1:6, n, replace = TRUE)
dat$guilt     <- clip(1.6 + 1.0 * dat$induction + 0.15 * (dat$eid - 5) + rnorm(n, 0, 0.8), 1, 5)
dat$pride     <- clip(2.8 - 0.8 * dat$induction + 0.15 * (dat$eid - 5) + rnorm(n, 0, 0.8), 1, 5)
ant <- dat$anticip == 1
dat$A_keep   <- ifelse(ant, clip(65 + 6 * dat$induction + rnorm(n, 0, 18), 0, 100), NA)
dat$A_donate <- ifelse(ant, clip(60 - 4 * dat$induction + rnorm(n, 0, 18), 0, 100), NA)
dat$M_keep   <- ifelse(ant, clip(round(5.0 + rnorm(n, 0, 1.2)), 1, 7), NA)
dat$M_donate <- ifelse(ant, clip(round(5.8 + 0.3 * dat$induction + rnorm(n, 0, 1.0)), 1, 7), NA)
lin <- qlogis(0.43) - 0.6 * dat$induction + 0.45 * (dat$guilt - 2) + 0.25 * (dat$eid - 5) +
       ifelse(ant, -0.025 * (dat$A_keep - dat$A_donate) + 0.3 * (dat$M_donate - dat$M_keep), 0)
dat$donate <- rbinom(n, 1, plogis(lin))

# =============================================================================
# PART 2. SCORES AND SAMPLES
# =============================================================================
z <- function(x) as.numeric(scale(x))
orig <- subset(dat, anticip == 0)                                   # H2-H4
antd <- subset(dat, anticip == 1)                                   # H5-H9

antd$dA <- z(z(antd$A_keep) - z(antd$A_donate))                     # keeping more pleasant than donating
antd$dM <- z(z(antd$M_donate) - z(antd$M_keep))                     # donating more right than keeping

dat$guilt_z  <- z(dat$guilt)
orig$guilt_z <- z(orig$guilt); orig$pride_z <- z(orig$pride)
for (v in c("eid", "age", "politics", "education")) {               # standardize covariates per sample
  dat[[paste0(v, "_z")]] <- z(dat[[v]]); orig[[paste0(v, "_z")]] <- z(orig[[v]]); antd[[paste0(v, "_z")]] <- z(antd[[v]])
}
cov_z <- "eid_z + gender + age_z + politics_z + education_z"

# =============================================================================
# PART 3. HELPERS
# =============================================================================
to_d <- function(b) oddsratio_to_d(b, log = TRUE)
lin_row <- function(H, claim, fit, term, dir) {
  ci <- confint(fit, term, level = 0.95)
  data.frame(H, claim, n = nobs(fit), dir, d = coef(fit)[[term]], lo = ci[1], hi = ci[2])
}
log_row <- function(H, claim, fit, term, dir) {
  ci <- confint.default(fit, term, level = 0.95)
  data.frame(H, claim, n = nobs(fit), dir, d = to_d(coef(fit)[[term]]), lo = to_d(ci[1]), hi = to_d(ci[2]))
}
outcome <- function(supported, not_supported) {
  if (supported && not_supported) "supported but too small to matter" else
  if (supported) "supported" else if (not_supported) "not supported" else "inconclusive"
}
decide <- function(r) {
  sup <- if (r$dir > 0) r$lo > 0 else r$hi < 0
  nsp <- if (r$dir > 0) r$hi < sesoi else r$lo > -sesoi
  r$outcome <- outcome(sup, nsp); r
}

# =============================================================================
# PART 4. THE FIVE MODELS
# =============================================================================

# Model 1: H1, all participants
m1 <- lm(as.formula(paste("guilt_z ~ induction +", cov_z)), data = dat)
r1 <- decide(lin_row("H1", "Being made to feel guilty raises felt guilt", m1, "induction", +1))

# Model 2: H2, original survey
m2 <- glm(as.formula(paste("donate ~ induction +", cov_z)), family = binomial, data = orig)
r2 <- decide(log_row("H2", "Being made to feel guilty lowers donation", m2, "induction", -1))

# Model 3: H3 and H4, original survey (PROCESS model 4, guilt and pride as parallel mediators)
p3 <- process(data = orig, y = "donate", x = "induction", m = c("guilt", "pride"), cov = covs,
              model = 4, boot = n_boot, seed = boot_seed, conf = 95, save = 1, progress = 0, outscreen = 0)
m3b <- glm(as.formula(paste("donate ~ induction + guilt_z + pride_z +", cov_z)), family = binomial, data = orig)
m3a <- lm(as.formula(paste("guilt_z ~ induction +", cov_z)), data = orig)
r3  <- decide(log_row("H3", "People who feel guiltier donate more", m3b, "guilt_z", +1))
link_a <- lin_row("H4 link a", "induction -> felt guilt", m3a, "induction", +1)
cov_raw <- paste(covs, collapse = " + ")
a_raw <- coef(lm(as.formula(paste("guilt ~ induction +", cov_raw)), data = orig))[["induction"]]
b_raw <- coef(glm(as.formula(paste("donate ~ induction + guilt + pride +", cov_raw)),
                  family = binomial, data = orig))[["guilt"]]
q <- quantile(p3$guilt_induction * p3$donate_guilt, c(0.025, 0.975))        # 95% percentile CI
r4 <- data.frame(H = "H4", claim = "Being made to feel guilty raises donation through felt guilt",
                 n = nobs(m3b), dir = +1, d = to_d(a_raw * b_raw), lo = to_d(q[[1]]), hi = to_d(q[[2]]))
r4$outcome <- outcome(r4$lo > 0, link_a$hi < sesoi || r3$hi < sesoi)

# Model 4: H5 and H6, anticipation measures
m5 <- lm(as.formula(paste("dA ~ induction +", cov_z)), data = antd)
m6 <- lm(as.formula(paste("dM ~ induction +", cov_z)), data = antd)
r5 <- decide(lin_row("H5", "Being made to feel guilty makes keeping more pleasant than donating", m5, "induction", +1))
r6 <- decide(lin_row("H6", "Being made to feel guilty makes donating more right than keeping", m6, "induction", +1))

# Model 5: H7, H8 and H9, anticipation measures
m7 <- glm(as.formula(paste("donate ~ induction + dA + dM +", cov_z)), family = binomial, data = antd)
r7 <- decide(log_row("H7", "Keeping more pleasant goes with less donation", m7, "dA", -1))
r8 <- decide(log_row("H8", "Donating more right goes with more donation", m7, "dM", +1))
V <- vcov(m7); s <- coef(m7)[["dA"]] + coef(m7)[["dM"]]
se_s <- sqrt(V["dA", "dA"] + V["dM", "dM"] + 2 * V["dA", "dM"])
r9 <- decide(data.frame(H = "H9", claim = "Pleasantness weighs more than rightness", n = nobs(m7), dir = -1,
                        d = to_d(s), lo = to_d(s - 1.96 * se_s), hi = to_d(s + 1.96 * se_s)))
if (!(r7$d < 0 && r8$d > 0)) r9$outcome <- paste(r9$outcome, "(not readable: H7 or H8 not in the predicted direction)")

# =============================================================================
# PART 5. RESULTS
# =============================================================================
res <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9)
out <- data.frame(H = res$H, Claim = res$claim, N = res$n, Prediction = ifelse(res$dir > 0, "+", "-"),
                  d = sprintf("%.2f", res$d), `95% CI` = sprintf("[%.2f, %.2f]", res$lo, res$hi),
                  Outcome = res$outcome, check.names = FALSE)
cat("\nPreregistered tests (Cohen's d, 95% CIs, smallest effect size of interest d = 0.20)\n\n")
print(out, row.names = FALSE, right = FALSE)
write.csv(out, file.path(root, "05_analysis/preregistered_analysis/analysis_results.csv"), row.names = FALSE)
cat("\nSaved 05_analysis/preregistered_analysis/analysis_results.csv\n")
sessionInfo()
