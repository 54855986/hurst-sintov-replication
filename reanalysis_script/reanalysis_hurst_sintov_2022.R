# =============================================================================
# Reanalysis of Hurst and Sintov (2022, Study 2) with the planned models
# =============================================================================
#
# PURPOSE
# Estimates H1 to H4 in the original data, with the same models, covariates and
# 95% confidence intervals as the planned analysis, so that the original effects
# can be compared with ours.
#
# DATA
# The raw data shared by the authors (Raw data.csv). We use the two conditions we
# replicate: control (Cond = 0) and experienced negative feedback (Cond = 1), and
# only participants the authors kept (Exclude = 0).
#
# MODELS (covariates: environmental self-identity, gender, age, political
# orientation and education, as in Hurst and Sintov, 2022, pp. 6, 8)
#   H1  felt guilt ~ induction                              linear regression
#   H2  donation ~ induction                                logistic regression
#   H3  donation ~ felt guilt (with pride and induction)    logistic regression
#   H4  induction -> felt guilt -> donation                 PROCESS model 4, guilt and pride
#                                                           as parallel mediators, 5,000 bootstrap samples
# Effects are Cohen's d. Continuous variables are standardized; log odds are
# converted with effectsize::oddsratio_to_d(log = TRUE).
#
# HOW TO RUN
#   PROCESS for R must be saved as 05_analysis/PROCESS/process.R (processmacro.org).
#   Install effectsize once: install.packages("effectsize"). Then run the whole script.
# =============================================================================

library(effectsize)

# ---- Paths ----------------------------------------------------------------------
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

# ---- Data ------------------------------------------------------------------------
raw <- read.csv(file.path(root, "05_analysis/reanalysis/Hurst_and_Sintov_2022/original_materials/Raw data.csv"),
                check.names = FALSE)
raw <- raw[raw$Exclude == 0 & raw$Cond %in% c(0, 1), ]
dat <- data.frame(induction = as.numeric(raw$Cond),       # 1 = experienced negative feedback, 0 = control
                  guilt     = as.numeric(raw$Guilt),
                  pride     = as.numeric(raw$Pride),
                  donate    = as.numeric(raw$Choice_2),   # 1 = donate, 0 = keep
                  eid       = as.numeric(raw$EI_Scale),
                  gender    = as.numeric(raw$Gend),
                  age       = as.numeric(raw$Age),
                  politics  = as.numeric(raw$Pol),
                  education = as.numeric(raw$Edu))
dat <- na.omit(dat)
cat("Participants with complete data:", nrow(dat), "\n")
cat(sprintf("Donated: control %.1f%%, guilt induction %.1f%%\n",
            100 * mean(dat$donate[dat$induction == 0]), 100 * mean(dat$donate[dat$induction == 1])))

# Standardized versions for the d scale
z <- function(x) as.numeric(scale(x))
for (v in c("guilt", "pride", "eid", "age", "politics", "education")) dat[[paste0(v, "_z")]] <- z(dat[[v]])
covs   <- c("eid", "gender", "age", "politics", "education")           # for PROCESS (raw units)
covs_z <- "eid_z + gender + age_z + politics_z + education_z"          # for the d-scale models

# ---- Helpers ---------------------------------------------------------------------
lin_row <- function(H, fit, term) {
  ci <- confint(fit, term, level = 0.95)
  data.frame(H, d = coef(fit)[[term]], lo = ci[1], hi = ci[2])
}
log_row <- function(H, fit, term) {
  ci <- confint.default(fit, term, level = 0.95)
  data.frame(H, d = oddsratio_to_d(coef(fit)[[term]], log = TRUE),
             lo = oddsratio_to_d(ci[1], log = TRUE), hi = oddsratio_to_d(ci[2], log = TRUE))
}

# ---- H1 to H3 ---------------------------------------------------------------------
m1 <- lm(as.formula(paste("guilt_z ~ induction +", covs_z)), data = dat)
m2 <- glm(as.formula(paste("donate ~ induction +", covs_z)), family = binomial, data = dat)
m3 <- glm(as.formula(paste("donate ~ induction + guilt_z + pride_z +", covs_z)), family = binomial, data = dat)
tab <- rbind(lin_row("H1 induction -> felt guilt", m1, "induction"),
             log_row("H2 induction -> donation", m2, "induction"),
             log_row("H3 felt guilt -> donation", m3, "guilt_z"))

# ---- H4: PROCESS model 4 ----------------------------------------------------------------
p4 <- process(data = dat, y = "donate", x = "induction", m = c("guilt", "pride"), cov = covs,
              model = 4, boot = 5000, seed = 2000000, conf = 95, save = 1, progress = 0, outscreen = 0)
a  <- coef(lm(as.formula(paste("guilt ~ induction +", paste(covs, collapse = " + "))), data = dat))[["induction"]]
b  <- coef(glm(as.formula(paste("donate ~ induction + guilt + pride +", paste(covs, collapse = " + "))),
              family = binomial, data = dat))[["guilt"]]
q  <- quantile(p4$guilt_induction * p4$donate_guilt, c(0.025, 0.975))
h4 <- data.frame(H = "H4 indirect effect through felt guilt (log odds)", d = a * b, lo = q[[1]], hi = q[[2]])

# ---- Results ----------------------------------------------------------------------------
out <- rbind(tab, h4)
out[, 2:4] <- round(out[, 2:4], 2)
cat("\nReanalysis of Hurst and Sintov (2022, Study 2), control vs. experienced negative feedback\n")
cat("H1-H3 as Cohen's d, H4 in log odds, all with 95% CIs\n\n")
print(out, row.names = FALSE)
