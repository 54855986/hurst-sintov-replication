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
#   rating     1 = rated keeping and donating $25 in a hypothetical version of the choice,
#                  after the intentions scale and before the raffle description (rating groups),
#              0 = followed the procedure of Hurst and Sintov (2022, Study 2) without these
#                  ratings (no-rating groups)
# The no-rating groups (rating = 0) carry H2 to H4; the rating groups (rating = 1) carry
# H5 to H9; H1 uses all participants.
#
# FIVE MODELS
#   Model 1  H1      all participants   linear:   experienced guilt ~ induction
#   Model 2  H2      no-rating groups   logistic: donation ~ induction
#   Model 3  H3, H4  no-rating groups   PROCESS model 4: induction -> experienced guilt and pride (parallel) -> donation
#   Model 4  H5, H6  rating groups      linear:   anticipated-affect difference ~ induction; moral-evaluation difference ~ induction
#   Model 5  H7-H9   rating groups      logistic: donation ~ induction + anticipated-affect difference + moral-evaluation difference
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
#   choice is more strongly related to the anticipated-affect difference than to the
#   moral-evaluation difference. It can be read
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
covs      <- c("eid", "woman", "gender_other", "age", "politics", "education")  # covariates (gender as two indicators, man = reference)
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
# exclude participants who failed the comprehension item and second responses of
# duplicate Prolific IDs (keep the first), and
# make sure these variables exist:
#   induction, rating    0/1
#   guilt, pride         mean of three items each (1-5)
#   A_keep, A_donate     how pleasant keeping / donating $25 would feel (0-100), rating = 1 only
#   M_keep, M_donate     how wrong or right keeping / donating $25 would be (1-7), rating = 1 only
#   donate               1 = donate, 0 = keep
#   eid                  environmental self-identity, mean of three items (1-7)
#   gender               1 = woman, 2 = man, 3 = another gender or not stated
#   age                  years
#   politics             political orientation (1 = very liberal to 7 = very conservative)
#   education            highest education (1-6)

n_norating <- 700; n_rating <- 500                                  # participants per group
dat <- data.frame(induction = c(rep(0:1, each = n_norating), rep(0:1, each = n_rating)),
                  rating    = rep(0:1, times = c(2 * n_norating, 2 * n_rating)))
n <- nrow(dat)
clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)
dat$eid       <- clip(round(rnorm(n, 5, 1.2), 1), 1, 7)
dat$gender    <- sample(1:3, n, replace = TRUE, prob = c(0.49, 0.49, 0.02))
dat$age       <- round(runif(n, 18, 75))
dat$politics  <- sample(1:7, n, replace = TRUE)
dat$education <- sample(1:6, n, replace = TRUE)
dat$guilt     <- clip(1.6 + 1.0 * dat$induction + 0.15 * (dat$eid - 5) + rnorm(n, 0, 0.8), 1, 5)
dat$pride     <- clip(2.8 - 0.8 * dat$induction + 0.15 * (dat$eid - 5) + rnorm(n, 0, 0.8), 1, 5)
ant <- dat$rating == 1
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
dat$woman        <- as.numeric(dat$gender == 1)                     # gender indicators, man = reference
dat$gender_other <- as.numeric(dat$gender == 3)

# Complete cases on the variables of each sample's models, before z-scoring
h1       <- dat[complete.cases(dat[, c("guilt", "induction", covs)]), ]                       # H1
norating <- subset(dat, rating == 0)                                                           # no-rating groups: H2-H4
norating <- norating[complete.cases(norating[, c("donate", "induction", "guilt", "pride", covs)]), ]
rated    <- subset(dat, rating == 1)                                                           # rating groups: H5-H9
rated    <- rated[complete.cases(rated[, c("donate", "induction", "A_keep", "A_donate", "M_keep", "M_donate", covs)]), ]

# Each person's own comparison, in the units of the measure (raw 0 = both options rated the same),
# then z-scored once for the models
rated$dA_raw <- rated$A_keep - rated$A_donate                       # slider points: keeping more pleasant than donating
rated$dM_raw <- rated$M_donate - rated$M_keep                       # scale points: donating more right than keeping
rated$dA <- z(rated$dA_raw)
rated$dM <- z(rated$dM_raw)

h1$guilt_z  <- z(h1$guilt)
norating$guilt_z <- z(norating$guilt); norating$pride_z <- z(norating$pride)
for (v in c("eid", "age", "politics", "education")) {               # standardize covariates per sample
  h1[[paste0(v, "_z")]] <- z(h1[[v]]); norating[[paste0(v, "_z")]] <- z(norating[[v]]); rated[[paste0(v, "_z")]] <- z(rated[[v]])
}
cov_z <- "eid_z + woman + gender_other + age_z + politics_z + education_z"

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
m1 <- lm(as.formula(paste("guilt_z ~ induction +", cov_z)), data = h1)
r1 <- decide(lin_row("H1", "Being made to feel guilty raises experienced guilt", m1, "induction", +1))

# Model 2: H2, no-rating groups
m2 <- glm(as.formula(paste("donate ~ induction +", cov_z)), family = binomial, data = norating)
if (!m2$converged) warning("Model 2 (H2) did not converge")
r2 <- decide(log_row("H2", "Being made to feel guilty lowers donation", m2, "induction", -1))

# Model 3: H3 and H4, no-rating groups (PROCESS model 4, guilt and pride as parallel mediators)
p3 <- process(data = norating, y = "donate", x = "induction", m = c("guilt", "pride"), cov = covs,
              model = 4, boot = n_boot, seed = boot_seed, conf = 95, save = 1, progress = 0, outscreen = 0)
m3b <- glm(as.formula(paste("donate ~ induction + guilt_z + pride_z +", cov_z)), family = binomial, data = norating)
if (!m3b$converged) warning("Model 3 (H3) did not converge")
m3a <- lm(as.formula(paste("guilt_z ~ induction +", cov_z)), data = norating)
r3  <- decide(log_row("H3", "People who feel guiltier donate more", m3b, "guilt_z", +1))
link_a <- lin_row("H4 link a", "induction -> experienced guilt", m3a, "induction", +1)
cov_raw <- paste(covs, collapse = " + ")
a_raw <- coef(lm(as.formula(paste("guilt ~ induction +", cov_raw)), data = norating))[["induction"]]
b_raw <- coef(glm(as.formula(paste("donate ~ induction + guilt + pride +", cov_raw)),
                  family = binomial, data = norating))[["guilt"]]
q <- quantile(p3$guilt_induction * p3$donate_guilt, c(0.025, 0.975))        # 95% percentile CI
r4 <- data.frame(H = "H4", claim = "The part of the effect of the guilt induction that runs through experienced guilt points toward donating",
                 n = nobs(m3b), dir = +1, d = to_d(a_raw * b_raw), lo = to_d(q[[1]]), hi = to_d(q[[2]]))
r4$outcome <- outcome(r4$lo > 0, link_a$hi < sesoi || r3$hi < sesoi)

# Model 4: H5 and H6, rating groups
m5 <- lm(as.formula(paste("dA ~ induction +", cov_z)), data = rated)
m6 <- lm(as.formula(paste("dM ~ induction +", cov_z)), data = rated)
r5 <- decide(lin_row("H5", "Being made to feel guilty raises the difference in anticipated affect (keeping minus donating)", m5, "induction", +1))
r6 <- decide(lin_row("H6", "Being made to feel guilty raises the difference in moral evaluation (donating minus keeping)", m6, "induction", +1))

# Model 5: H7, H8 and H9, rating groups
m7 <- glm(as.formula(paste("donate ~ induction + dA + dM +", cov_z)), family = binomial, data = rated)
if (!m7$converged) warning("Model 5 (H7-H9) did not converge")
r7 <- decide(log_row("H7", "A larger difference in anticipated affect goes with less donation", m7, "dA", -1))
r8 <- decide(log_row("H8", "A larger difference in moral evaluation goes with more donation", m7, "dM", +1))
V <- vcov(m7); s <- coef(m7)[["dA"]] + coef(m7)[["dM"]]
se_s <- sqrt(V["dA", "dA"] + V["dM", "dM"] + 2 * V["dA", "dM"])
r9 <- decide(data.frame(H = "H9", claim = "The choice is more strongly related to the difference in anticipated affect than to the difference in moral evaluation", n = nobs(m7), dir = -1,
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
# Descriptives of the ratings (not tests). "Equal" means a raw difference of exactly 0.
grp <- factor(rated$induction, levels = 0:1, labels = c("control", "guilt"))
cat("\nRatings in the rating groups, by guilt induction (descriptives, not tests)\n")
for (g in levels(grp)) {
  s <- rated[grp == g, ]
  cat("\n", g, "\n", sep = "")
  cat(sprintf("  A_keep    M = %.2f, SD = %.2f\n", mean(s$A_keep, na.rm = TRUE), sd(s$A_keep, na.rm = TRUE)))
  cat(sprintf("  A_donate  M = %.2f, SD = %.2f\n", mean(s$A_donate, na.rm = TRUE), sd(s$A_donate, na.rm = TRUE)))
  cat(sprintf("  M_keep    M = %.2f, SD = %.2f\n", mean(s$M_keep, na.rm = TRUE), sd(s$M_keep, na.rm = TRUE)))
  cat(sprintf("  M_donate  M = %.2f, SD = %.2f, at 7 = %.0f%%\n", mean(s$M_donate, na.rm = TRUE),
              sd(s$M_donate, na.rm = TRUE), 100 * mean(s$M_donate == 7, na.rm = TRUE)))
  cat(sprintf("  A_keep - A_donate  M = %.2f, SD = %.2f, below 0 = %.0f%%, at 0 = %.0f%%, above 0 = %.0f%%\n",
              mean(s$dA_raw, na.rm = TRUE), sd(s$dA_raw, na.rm = TRUE), 100 * mean(s$dA_raw < 0, na.rm = TRUE),
              100 * mean(s$dA_raw == 0, na.rm = TRUE), 100 * mean(s$dA_raw > 0, na.rm = TRUE)))
  cat(sprintf("  M_donate - M_keep  M = %.2f, SD = %.2f, below 0 = %.0f%%, at 0 = %.0f%%, above 0 = %.0f%%\n",
              mean(s$dM_raw, na.rm = TRUE), sd(s$dM_raw, na.rm = TRUE), 100 * mean(s$dM_raw < 0, na.rm = TRUE),
              100 * mean(s$dM_raw == 0, na.rm = TRUE), 100 * mean(s$dM_raw > 0, na.rm = TRUE)))
}

sessionInfo()
