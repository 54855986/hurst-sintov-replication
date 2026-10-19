# Power Analysis for H1 to H8

2026-10-19

## About this document

The sample size of the study was set by a fixed budget, which is a
justification by resource constraints (Lakens, 2022). The budget buys
about 2,800 participants, about 2,660 after exclusions. This document
checks by simulation that this sample is informative for every
confirmatory hypothesis. For each hypothesis it estimates two
probabilities.

- **Power to find** is the chance that the test detects the effect when
  the effect exists.
- **Power to show absent** is the chance that the test shows the effect
  to be absent, or too small to matter, when it does not exist.

The document also shows how the split of participants between the two
parts of the study was chosen, and which effects the design can detect
(sensitivity analysis).

**Design.**

| Cell | Feedback | Asked about anticipated emotions | Part | Participants after exclusions |
|----|----|----|----|----|
| 1 | None | No | Replication | 798 |
| 2 | Negative | No | Replication | 798 |
| 3 | None | Yes | Extension | 532 |
| 4 | Negative | Yes | Extension | 532 |

## Settings

``` r
n_not_asked  <- 798               # participants per cell in the replication part
n_asked      <- 532               # participants per cell in the extension part
n_sim        <- 5000              # simulated data sets per world and donation rate (final: 5000)
n_sim_scan   <- 5000               # simulated data sets per allocation (final: 5000)
p_control    <- c(.36, .43, .50)  # donation rate in the control group
d_donation   <- -0.19             # total effect of feedback on donation (H2)
d_to_logodds <- pi / sqrt(3)      # converts d to log-odds
z90          <- qnorm(.95)        # critical value for a 90% CI
# Results are saved and reused, so the document renders quickly after the first run.
# Delete the file, or change a setting above, to run the simulation again.
cache_file <- file.path(project_root, "05_analysis/power_analysis/power_results_simulated.rds")
```

- **Donation rate.** In the control condition of the original data, 43%
  donated (83 of 193). The simulation uses this rate and the limits of
  its 95% confidence interval, 36% and 50%. Power is lowest when the
  rate lies furthest from 50%.
- **Expected drop in donation.** d = −0.19, smaller than the estimate
  from our reanalysis of the original data (d = −0.27), because a single
  estimate found after the fact is likely to be too large.
- **Number of simulated data sets.** 5000 per scenario for the main
  results and 5000 per allocation. With 1,000 data sets, a power value
  is precise to about ±.02.

## Hypotheses and decision rules

| H | Claim | Cells | Limit |
|----|----|----|----|
| H1 | Negative feedback raises guilt | 1–4 | d = 0.155 |
| H2 | Negative feedback lowers donation | 1–2 | d = 0.189 |
| H3 | Path from feedback through guilt to donation | 1–2 | d = 0.20 per link |
| H4 | Negative direct effect of feedback on donation | 1–2 | d = 0.20 |
| H5 | No effect of feedback on intentions (equivalence) | 1–4 | ±0.155 |
| H6 | Feedback makes keeping feel better | 3–4 | d = 0.20 |
| H7 | Keeping feels better, less donation | 3–4 | d = 0.20 per SD |
| H8 | Paths from feedback through anticipated emotions to donation (serial model: guilt and pride, then anticipated emotions) | 3–4 | links H6, H7 |

The limits 0.155 and 0.189 are the effects the original study had a 33%
chance of detecting. The limit 0.20 is the smallest effect of interest
for the other tests. A test **finds** an effect when the 90% confidence
interval lies on the predicted side of zero, and **shows it absent**
when the 90% interval lies on the near side of the limit. H5 is shown
absent when its 90% interval lies inside ±0.155. H3 counts as found when
both of its links are found. H8 is the sum of all paths through
anticipated emotions in a serial model (PROCESS model 6 in the analysis)
and counts as found when its 90% Monte Carlo interval lies below zero.
Both paths count as absent when one of their links is absent.

## The simulated worlds

Every scenario is simulated in three worlds.

- **Effects.** All effects exist. Feedback raises guilt (d = 1.31) and
  lowers pride, lowers donation in total (d = −0.19) and shifts
  anticipated emotions toward keeping (d = about 0.20), partly directly
  and partly because the guilt it evokes changes anticipated emotions.
  Guilt raises donation, and anticipated emotions toward keeping lower
  it (d = −0.20 per SD). Feedback has no effect on intentions. This
  world gives the **power to find**.
- **Null.** Feedback still raises guilt and lowers pride, as in every
  earlier study, but has no other effect, and neither guilt nor
  anticipated emotions are related to donation. This world gives the
  **power to show absent** for H2 to H8. Keeping the guilt effect
  matters for H4, because guilt and feedback overlap in the model, which
  makes the direct effect harder to pin down.
- **No guilt.** As the null world, and feedback does not change guilt
  either. This world gives the power to show that H1 is absent.

``` r
# One simulated data set. n_cells: participants in cells 1-4; b0: intercept of donation
# (log-odds); cprime: direct effect of feedback on donation (log-odds).
simulate <- function(n_cells, world, b0, cprime) {
  eff <- world == "effects"
  emo <- world != "no_guilt"
  X <- rep(c(0, 1, 0, 1), n_cells)                         # feedback
  S <- rep(c(0, 0, 1, 1), n_cells)                         # asked about anticipated emotions
  N <- length(X)
  E <- rnorm(N)                                            # environmental self-identity
  e <- mvrnorm(N, c(0, 0), matrix(c(1, -.3, -.3, 1), 2))
  guilt <- 1.31 * X * emo + .35 * E + e[, 1]
  pride <- -1.29 * X * emo - .20 * E + e[, 2]
  pei   <- .30 * E + rnorm(N, 0, sqrt(1 - .30^2))          # intentions, no effect of feedback
  A     <- (.095 * X + .10 * guilt / sd(guilt)) * eff + .20 * E + rnorm(N)              # anticipated-emotions score
  b_guilt <- if (eff) .486 else 0                          # log-odds per SD of guilt
  b_A     <- if (eff) -.20 * d_to_logodds else 0           # log-odds per SD of A
  lin <- b0 + cprime * X + .725 * E + b_guilt * guilt / sd(guilt) + b_A * A / sd(A)
  data.frame(X, S, E, guilt, pride, pei, A, donate = rbinom(N, 1, plogis(lin)))
}

# Calibration on one large data set: b0 gives the control group the donation rate p, and
# (effects world) cprime makes the total effect of feedback, controlling for E, equal d_donation.
calibrate <- function(world, p) {
  big <- rep(25000, 4)
  rate_gap <- function(b0, cprime) {
    set.seed(1); d <- simulate(big, world, b0, cprime); mean(d$donate[d$X == 0]) - p
  }
  total_gap <- function(cprime) {
    b0 <- uniroot(rate_gap, c(-5, 5), cprime = cprime)$root
    set.seed(1); d <- simulate(big, world, b0, cprime)
    coef(glm(donate ~ X + E, binomial, d))["X"] / d_to_logodds - d_donation
  }
  cprime <- if (world == "effects") uniroot(total_gap, c(-3, 1))$root else 0
  c(b0 = uniroot(rate_gap, c(-5, 5), cprime = cprime)$root, cprime = cprime)
}

# Tests on one data set; each helper returns the 90% CI of an effect in d units.
ci90  <- function(est, se) unname(c(est - z90 * se, est + z90 * se))
lin_d <- function(fit, term) { s <- summary(fit)$coefficients[term, 1:2]; ci90(s[1], s[2]) }
log_d <- function(fit, term) { s <- summary(fit)$coefficients[term, 1:2] / d_to_logodds; ci90(s[1], s[2]) }
# Path through the anticipated-emotions score in the serial model (guilt, pride, score),
# with a Monte Carlo 90% CI (2,000 draws) as a fast stand-in for the PROCESS bootstrap.
serial_ci90 <- function(d34) {
  z <- function(v) as.numeric(scale(v))
  d34$gz <- z(d34$guilt); d34$pz <- z(d34$pride); d34$az <- z(d34$A)
  draw <- function(fit, terms) mvrnorm(2000, coef(fit)[terms], vcov(fit)[terms, terms, drop = FALSE])
  g <- rnorm(2000, coef(mg <- lm(gz ~ X + E, d34))["X"], sqrt(vcov(mg)["X", "X"]))
  p <- draw(lm(pz ~ X + gz + E, d34), c("X", "gz"))
  a <- draw(lm(az ~ X + gz + pz + E, d34), c("X", "gz", "pz"))
  my <- glm(donate ~ X + gz + pz + az + E, binomial, d34)
  b <- rnorm(2000, coef(my)["az"], sqrt(vcov(my)["az", "az"]))
  ind <- b * (a[, "X"] + g * a[, "gz"] + p[, "X"] * a[, "pz"] + g * p[, "gz"] * a[, "pz"])
  unname(quantile(ind / d_to_logodds, c(.05, .95)))
}

run_tests <- function(d) {
  d12 <- d[d$S == 0, ]; d34 <- d[d$S == 1, ]
  z <- function(v) as.numeric(scale(v))
  h1 <- lin_d(lm(z(guilt) ~ X + E, d), "X")
  h2 <- log_d(glm(donate ~ X + E, binomial, d12), "X")
  a3 <- lin_d(lm(z(guilt) ~ X + E, d12), "X")
  m3 <- glm(donate ~ X + z(guilt) + z(pride) + E, binomial, d12)
  b3 <- log_d(m3, "z(guilt)"); h4 <- log_d(m3, "X")
  h5 <- lin_d(lm(z(pei) ~ X + E, d), "X")
  h6 <- lin_d(lm(z(A) ~ X + E, d34), "X")
  h7 <- log_d(glm(donate ~ z(A) + X + z(guilt) + z(pride) + E, binomial, d34), "z(A)")
  c(H1_found = h1[1] > 0,  H1_absent = h1[2] < .155,
    H2_found = h2[2] < 0,  H2_absent = h2[1] > -.189,
    H3_found = a3[1] > 0 & b3[1] > 0, H3_absent = a3[2] < .20 | b3[2] < .20,
    H4_found = h4[2] < 0,  H4_absent = h4[1] > -.20,
    H5_absent = h5[1] > -.155 & h5[2] < .155,
    H6_found = h6[1] > 0,  H6_absent = h6[2] < .20,
    H7_found = h7[2] < 0,  H7_absent = h7[1] > -.20,
    H8_found = serial_ci90(d34)[2] < 0, H8_absent = h6[2] < .20 | h7[1] > -.20)
}
```

``` r
settings <- list(n_not_asked, n_asked, n_sim, n_sim_scan, p_control, d_donation, "serial H8")
if (file.exists(cache_file) && identical(readRDS(cache_file)$settings, settings)) {
  saved <- readRDS(cache_file); results <- saved$results; scan <- saved$scan
} else {
  worlds <- c("effects", "null", "no_guilt")
  cals <- list(); results <- list(); run <- 0
  n_cells <- c(n_not_asked, n_not_asked, n_asked, n_asked)
  for (p in p_control) for (w in worlds) {
    cals[[paste(w, p)]] <- calibrate(w, p)
    run <- run + 1; set.seed(2026 + run)                   # calibrate() resets the seed
    cal <- cals[[paste(w, p)]]
    results[[paste(w, p)]] <- rowMeans(replicate(n_sim, run_tests(simulate(n_cells, w, cal["b0"], cal["cprime"]))))
  }
  # Allocations: the total is fixed, the split between the two parts varies.
  ratios <- list("1:1" = c(1, 1), "5:4" = c(5, 4), "4:3" = c(4, 3), "3:2" = c(3, 2), "2:1" = c(2, 1))
  per_feedback <- n_not_asked + n_asked
  scan <- do.call(rbind, lapply(names(ratios), function(r) {
    w12 <- ratios[[r]]; n12 <- round(per_feedback * w12[1] / sum(w12)); n34 <- per_feedback - n12
    do.call(rbind, lapply(p_control, function(p) {
      sim <- function(w) { cal <- cals[[paste(w, p)]]; set.seed(3000 + n12)
        rowMeans(replicate(n_sim_scan, run_tests(simulate(c(n12, n12, n34, n34), w, cal["b0"], cal["cprime"])))) }
      e <- sim("effects"); nl <- sim("null")
      data.frame(ratio = r, n12 = n12, n34 = n34, p = p, H2_find = e[["H2_found"]],
                 H4_absent = nl[["H4_absent"]], H6_find = e[["H6_found"]], H8_find = e[["H8_found"]])
    }))
  }))
  saveRDS(list(settings = settings, results = results, scan = scan), cache_file)
}
```

## Power for each hypothesis

``` r
hyp <- paste0("H", 1:8)
long <- do.call(rbind, lapply(p_control, function(p) {
  data.frame(H = hyp, p = p,
             find   = sapply(hyp, function(h) results[[paste("effects", p)]][paste0(h, "_found")]),
             absent = sapply(hyp, function(h) {
               w <- if (h == "H1") "no_guilt" else "null"; results[[paste(w, p)]][paste0(h, "_absent")] }))
}))
rng <- function(v) if (all(is.na(v))) "-" else {
  r <- round(range(v, na.rm = TRUE), 2); if (r[1] == r[2]) sprintf("%.2f", r[1]) else sprintf("%.2f to %.2f", r[1], r[2])
}
tab <- data.frame(Hypothesis = hyp,
                  Cells = c("1-4", "1-2", "1-2", "1-2", "1-4", "3-4", "3-4", "3-4"),
                  `Power to find` = sapply(hyp, function(h) rng(long$find[long$H == h])),
                  `Power to show absent` = sapply(hyp, function(h) rng(long$absent[long$H == h])),
                  check.names = FALSE)
kable(tab, row.names = FALSE, caption = "Range over the three donation rates in the control group (36%, 43%, 50%)")
```

| Hypothesis | Cells | Power to find | Power to show absent |
|:-----------|:------|:--------------|:---------------------|
| H1         | 1-4   | 1.00          | 1.00                 |
| H2         | 1-2   | 0.92 to 0.95  | 0.92 to 0.94         |
| H3         | 1-2   | 1.00          | 1.00                 |
| H4         | 1-2   | 1.00          | 0.82 to 0.84         |
| H5         | 1-4   | \-            | 0.99                 |
| H6         | 3-4   | 0.94 to 0.95  | 0.95 to 0.96         |
| H7         | 3-4   | 1.00          | 1.00                 |
| H8         | 3-4   | 0.94 to 0.95  | 1.00                 |

Range over the three donation rates in the control group (36%, 43%, 50%)

``` r
autumn <- c(rust = "#A63D1E", orange = "#D9822B", gold = "#E3A33B",
            olive = "#6B7F3A", brown = "#5C3A21", plum = "#7A3B4F")
theme_set(theme_minimal(base_size = 12) +
          theme(plot.title = element_text(face = "bold", colour = "black"), axis.text = element_text(colour = "black"),
                legend.position = "bottom", panel.grid.minor = element_blank()))
pl <- rbind(transform(long[, c("H", "p")], kind = "Power to find", power = long$find),
            transform(long[, c("H", "p")], kind = "Power to show absent", power = long$absent))
pl <- pl[!is.na(pl$power), ]
pl$H <- factor(pl$H, levels = rev(hyp))
sumpl <- do.call(rbind, lapply(split(pl, list(pl$H, pl$kind), drop = TRUE), function(x)
  data.frame(H = x$H[1], kind = x$kind[1], lo = min(x$power), hi = max(x$power), mid = x$power[x$p == .43])))
ggplot(sumpl, aes(y = H, colour = kind)) +
  annotate("rect", xmin = .80, xmax = .90, ymin = -Inf, ymax = Inf, fill = autumn[["gold"]], alpha = .12) +
  geom_vline(xintercept = c(.80, .90), linetype = c("dashed", "solid"), colour = autumn[["olive"]]) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .3, linewidth = 1,
                position = position_dodge(.5)) +
  geom_point(aes(x = mid), size = 3, position = position_dodge(.5)) +
  scale_colour_manual(values = c("Power to find" = autumn[["rust"]], "Power to show absent" = autumn[["orange"]]),
                      name = NULL) +
  coord_cartesian(xlim = c(.6, 1)) +
  labs(x = "Power (point: 43% donation rate; bar: range from 36% to 50%)", y = NULL,
       title = "Power for H1 to H8",
       subtitle = "Dashed line: minimum of .80. Solid line: target of .90")
```

![](power_fixed_budget_simple_files/figure-commonmark/power-plot-1.png)

Every confirmatory hypothesis reaches the target of .90 for both kinds
of power, except the power to show that the negative direct effect (H4)
is absent, which stays above the minimum of .80. H5 predicts no effect,
so it has no power to find.

## Choice of the allocation

The total sample is fixed. H2 to H4 are tested only in the replication
part and H6 to H8 only in the extension part, so every participant added
to one part is taken from the other. Five ratios were compared. The
chosen ratio is the most even one at which the drop in donation (H2) and
the tests of anticipated emotions (H6, H8) reach .90 and the direct
effect (H4) can be shown absent with at least .80.

``` r
scan_min <- aggregate(cbind(H2_find, H4_absent, H6_find, H8_find) ~ ratio + n12 + n34, data = scan, FUN = min)
scan_min <- scan_min[match(c("1:1", "5:4", "4:3", "3:2", "2:1"), scan_min$ratio), ]
kable(transform(scan_min, H2_find = round(H2_find, 2), H4_absent = round(H4_absent, 2),
                H6_find = round(H6_find, 2), H8_find = round(H8_find, 2)),
      row.names = FALSE, col.names = c("Ratio", "Per cell, replication", "Per cell, extension",
                                       "H2 find", "H4 absent", "H6 find", "H8 find"),
      caption = "Lowest power over the three donation rates")
```

| Ratio | Per cell, replication | Per cell, extension | H2 find | H4 absent | H6 find | H8 find |
|:---|---:|---:|---:|---:|---:|---:|
| 1:1 | 665 | 665 | 0.87 | 0.75 | 0.98 | 0.98 |
| 5:4 | 739 | 591 | 0.90 | 0.80 | 0.96 | 0.96 |
| 4:3 | 760 | 570 | 0.90 | 0.81 | 0.96 | 0.96 |
| 3:2 | 798 | 532 | 0.92 | 0.82 | 0.95 | 0.95 |
| 2:1 | 887 | 443 | 0.94 | 0.86 | 0.90 | 0.90 |

Lowest power over the three donation rates

``` r
sl <- do.call(rbind, lapply(c("H2_find", "H4_absent", "H6_find", "H8_find"), function(v)
  data.frame(ratio = scan_min$ratio, n12 = scan_min$n12, test = v, power = scan_min[[v]])))
sl$test <- factor(sl$test, levels = c("H2_find", "H4_absent", "H6_find", "H8_find"),
                  labels = c("H2 find (target .90)", "H4 show absent (target .80)",
                             "H6 find (target .90)", "H8 find (target .90)"))
chosen <- scan_min$n12[scan_min$ratio == "3:2"]
ggplot(sl, aes(n12, power, colour = test)) +
  annotate("rect", xmin = chosen - 8, xmax = chosen + 8, ymin = -Inf, ymax = Inf, fill = autumn[["gold"]], alpha = .25) +
  geom_hline(yintercept = c(.80, .90), linetype = c("dashed", "solid"), colour = autumn[["olive"]]) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
  scale_x_continuous(breaks = scan_min$n12, labels = paste0(scan_min$ratio, "\n", scan_min$n12, " / ", scan_min$n34)) +
  scale_colour_manual(values = unname(autumn[c("rust", "plum", "orange", "brown")]), name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Ratio and participants per cell (replication / extension)", y = "Lowest power",
       title = "Power as participants move from the extension to the replication part",
       subtitle = "Shaded: the chosen ratio 3:2")
```

![](power_fixed_budget_simple_files/figure-commonmark/allocation-plot-1.png)

Moving participants into the replication part raises the power for H2
and H4 and lowers it for H6 and H8. An equal split leaves H2 and H4
below their targets, and 3:2 is the first ratio at which all four tests
reach them.

## Sensitivity analysis

Following Lakens (2022), the last figure shows which effects the design
can detect. It uses a normal approximation with one-sided tests at an
alpha level of .05 and ignores the gain from the covariate, so the real
power is slightly higher.

``` r
se_lin <- function(n) sqrt(2 / n)                                     # n per group
se_log <- function(n, p) sqrt(2 / (n * p * (1 - p))) / d_to_logodds   # d units, donation
parts <- list("H1, H5 (all participants)" = se_lin(n_not_asked + n_asked),
              "H2 (replication part, donation 43%)" = se_log(n_not_asked, .43),
              "H6 (extension part)" = se_lin(n_asked))
sens <- t(sapply(parts, function(se) round((z90 + qnorm(c(.80, .90, .95))) * se, 2)))
kable(sens, col.names = c("Power .80", "Power .90", "Power .95"),
      caption = "Smallest effect (d) detected with each level of power")
```

|                                     | Power .80 | Power .90 | Power .95 |
|:------------------------------------|----------:|----------:|----------:|
| H1, H5 (all participants)           |      0.10 |      0.11 |      0.13 |
| H2 (replication part, donation 43%) |      0.14 |      0.16 |      0.18 |
| H6 (extension part)                 |      0.15 |      0.18 |      0.20 |

Smallest effect (d) detected with each level of power

``` r
dgrid <- seq(0, .4, .005)
curves <- do.call(rbind, lapply(names(parts), function(nm)
  data.frame(test = nm, d = dgrid, power = pnorm(dgrid / parts[[nm]] - z90))))
ggplot(curves, aes(d, power, colour = test)) +
  geom_hline(yintercept = c(.80, .90, .95), linetype = c("dashed", "solid", "dotted"), colour = autumn[["olive"]]) +
  geom_vline(xintercept = c(.155, .19, .20), colour = "grey70", linetype = "dotted") +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = unname(autumn[c("olive", "rust", "orange")]), name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "True effect (Cohen's d)", y = "Power to find the effect",
       title = "Sensitivity: power across true effect sizes",
       subtitle = "Horizontal lines: power .80, .90 and .95. Vertical lines: d = 0.155, 0.19 and 0.20")
```

![](power_fixed_budget_simple_files/figure-commonmark/sensitivity-1.png)

## Simplifications

Environmental self-identity is the only covariate in the simulation. The
other covariates of the preregistered models add little precision in a
randomized design, so the real power is slightly higher. H3 is tested by
the joint significance of its two links and H8 by a Monte Carlo
confidence interval (2,000 draws), both instead of the PROCESS
bootstrap, and confidence intervals are Wald intervals. Effects from
logistic models are converted to Cohen’s *d* with d = log-odds × √3/π.

``` r
sessionInfo()
```

    R version 4.5.2 (2025-10-31)
    Platform: aarch64-apple-darwin20
    Running under: macOS Sequoia 15.3.2

    Matrix products: default
    BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

    locale:
    [1] C

    time zone: Europe/Prague
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] knitr_1.51    ggplot2_4.0.3 MASS_7.3-65  

    loaded via a namespace (and not attached):
     [1] vctrs_0.7.1        cli_3.6.6          rlang_1.2.0        xfun_0.56         
     [5] otel_0.2.0         generics_0.1.4     S7_0.2.1           jsonlite_2.0.0    
     [9] labeling_0.4.3     glue_1.8.0         htmltools_0.5.9    scales_1.4.0      
    [13] rmarkdown_2.30     grid_4.5.2         tibble_3.3.1       evaluate_1.0.5    
    [17] fastmap_1.2.0      yaml_2.3.12        lifecycle_1.0.5    compiler_4.5.2    
    [21] dplyr_1.2.1        RColorBrewer_1.1-3 pkgconfig_2.0.3    rstudioapi_0.18.0 
    [25] farver_2.1.2       digest_0.6.39      R6_2.6.1           tidyselect_1.2.1  
    [29] pillar_1.11.1      magrittr_2.0.4     withr_3.0.2        tools_4.5.2       
    [33] gtable_0.3.6      
