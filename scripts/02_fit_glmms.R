#!/usr/bin/env Rscript
# ============================================================================
# Poisson and negative-binomial GLMMs with the same mean model as fcDPGLM
# (same design matrix X, exposure offset and rat random intercept), compared
# with the fcDPGLM fit: stage-specific mean curves m_d(t), the stage
# difference Delta(t), stage contrasts, and posterior predictive checks of the
# zero proportion, variance, upper percentiles and maximum.
#
# Usage (from the repository root, after scripts/01_fit_fcdpglm.R):
#   Rscript scripts/02_fit_glmms.R
# Optional arguments: dpglm_path=results/fcdpglm_fit.rds out_dir=results/glmm_comparison
# ============================================================================

suppressPackageStartupMessages({
  library(glmmTMB)
  library(splines)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(MASS)
  library(Rcpp)
})

set.seed(20260923)

## ===========================================================================
## SETTINGS
## ===========================================================================

dpglm_path <- "results/fcdpglm_fit.rds"
sampler_path <- "R/fcdpglm.R"

# Full-data uncertainty / predictive simulation
n_beta_sim <- 5000L
n_pp_sim   <- 1000L

# Output directory
out_dir <- "results/glmm_comparison"
# Optional command-line overrides: dpglm_path=..., out_dir=...
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) next
  if (kv[1] == "dpglm_path") dpglm_path <- kv[2]
  if (kv[1] == "out_dir") out_dir <- kv[2]
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## ===========================================================================
## HELPERS
## ===========================================================================

summ <- function(x) {
  q <- unname(quantile(x, c(.025, .5, .975), na.rm = TRUE))
  c(mean = mean(x, na.rm = TRUE),
    median = q[2],
    lower95 = q[1],
    upper95 = q[3])
}

rmvn_safe <- function(n, mu, Sigma) {
  Sigma <- (Sigma + t(Sigma)) / 2
  ee <- eigen(Sigma, symmetric = TRUE)
  ee$values <- pmax(ee$values, 1e-10)
  Sigma_pd <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
  MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma_pd)
}

make_glmm_df <- function(y, exposure, rat, X) {
  design_df <- as.data.frame(X)
  names(design_df) <- paste0("X_", make.names(colnames(X)))
  data.frame(
    y = y,
    exposure = exposure,
    rat = factor(rat),
    design_df,
    check.names = FALSE
  )
}

make_glmm_formula <- function(X) {
  xnames <- paste0("X_", make.names(colnames(X)))
  as.formula(
    paste(
      "y ~ 0 +", paste(xnames, collapse = " + "),
      "+ offset(log(exposure)) + (1 | rat)"
    )
  )
}

check_coef_alignment <- function(fit, X) {
  expected <- paste0("X_", make.names(colnames(X)))
  got <- names(fixef(fit)$cond)
  if (!identical(got, expected)) {
    stop("Coefficient/design mismatch.\nExpected: ",
         paste(expected, collapse = ", "),
         "\nGot: ", paste(got, collapse = ", "), call. = FALSE)
  }
}

# Summary statistics for the predictive checks.
pred_stats <- function(y) {
  c(
    mean = mean(y),
    variance = var(y),
    zero_rate = mean(y == 0),
    p90 = unname(quantile(y, .90)),
    p95 = unname(quantile(y, .95)),
    maximum = max(y)
  )
}

pp_summary <- function(y_obs, yrep) {
  obs <- pred_stats(y_obs)
  reps <- t(apply(yrep, 1L, pred_stats))
  data.frame(
    statistic = names(obs),
    observed = as.numeric(obs),
    predictive_median = apply(reps, 2, median),
    lower95 = apply(reps, 2, quantile, .025),
    upper95 = apply(reps, 2, quantile, .975),
    predictive_p = colMeans(sweep(reps, 2, obs, ">=")),
    row.names = NULL
  )
}

# Draw one categorical count per row from row-wise probabilities.
draw_counts <- function(prob, values) {
  cs <- t(apply(prob, 1L, cumsum))
  cs[, ncol(cs)] <- 1
  values[max.col(cs >= runif(nrow(prob)), ties.method = "first")]
}

tilted_pmf <- function(theta, weights, support) {
  lp <- outer(theta, support) + rep(log(weights), each = length(theta))
  lp <- lp - apply(lp, 1L, max)
  p <- exp(lp)
  p / rowSums(p)
}


## ===========================================================================
## COMPARISON
## ===========================================================================

obj <- readRDS(dpglm_path)
dat <- as.data.frame(obj$data)
X   <- as.matrix(obj$design_matrix)
fit_dp <- obj$rat_dpglm

stopifnot(
  nrow(dat) == nrow(X),
  length(dat$y) == nrow(X),
  all(dat$exposure > 0)
)

# X = [intercept, stage2, stage1 spline block, stage2 spline block], taken
# from the fcDPGLM fit so that both GLMMs use the same mean model.
glmm_dat <- make_glmm_df(dat$y, dat$exposure, dat$int_id, X)
glmm_formula <- make_glmm_formula(X)

cat("\nFitting full-data Poisson GLMM...\n")
fit_pois <- glmmTMB(
  glmm_formula,
  data = glmm_dat,
  family = poisson(link = "log")
)

cat("Fitting full-data NB2 GLMM...\n")
fit_nb <- glmmTMB(
  glmm_formula,
  data = glmm_dat,
  family = nbinom2(link = "log")
)

check_coef_alignment(fit_pois, X)
check_coef_alignment(fit_nb, X)

saveRDS(
  list(poisson = fit_pois, negative_binomial = fit_nb),
  file.path(out_dir, "full_data_glmms.rds")
)


## ===========================================================================
## A1. COMMON GRID FOR m_d(t) AND Delta(t)
## ===========================================================================

cn <- colnames(X)
idx_int   <- which(cn == "intercept")
idx_stage <- which(cn == "stage2")
idx_g1    <- grep("^stage1_spline", cn)
idx_g2    <- grep("^stage2_spline", cn)

stopifnot(length(idx_int) == 1L, length(idx_stage) == 1L,
          length(idx_g1) == length(idx_g2), length(idx_g1) > 0L)

# Spline basis as in scripts/01_fit_fcdpglm.R.
basis <- ns(dat$t, df = 3, intercept = FALSE)
X_check <- cbind(
  intercept = 1,
  stage2 = as.integer(dat$stage == 2L),
  predict(basis, dat$t) * (dat$stage == 1L),
  predict(basis, dat$t) * (dat$stage == 2L)
)
colnames(X_check) <- cn
stopifnot(max(abs(X_check - X)) < 1e-8)

t_of_day <- function(day, stage_value) {
  dd <- dat$day[dat$stage == stage_value]
  (day - mean(dd)) / sd(dd)
}

rng <- sapply(split(dat$day, dat$stage), range)
day_grid <- seq(max(rng[1, ]), min(rng[2, ]))
G <- length(day_grid)

B1 <- predict(basis, t_of_day(day_grid, 1L))
B2 <- predict(basis, t_of_day(day_grid, 2L))

Xg1 <- matrix(0, G, ncol(X), dimnames = list(NULL, cn))
Xg2 <- matrix(0, G, ncol(X), dimnames = list(NULL, cn))
Xg1[, idx_int] <- 1
Xg2[, idx_int] <- 1
Xg2[, idx_stage] <- 1
Xg1[, idx_g1] <- B1
Xg2[, idx_g2] <- B2

# Typical rat (a = 0) with one day of exposure, as in the definition of
# m_d(t).


## ===========================================================================
## A2. DPGLM CURVES
## ===========================================================================

burn <- if (!is.null(fit_dp$tuning$adapt_until)) {
  as.integer(fit_dp$tuning$adapt_until)
} else {
  floor(nrow(fit_dp$beta) / 2)
}
keep <- seq.int(burn + 1L, nrow(fit_dp$beta))
beta_dp <- fit_dp$beta[keep, , drop = FALSE]
colnames(beta_dp) <- cn

eta1_dp <- beta_dp %*% t(Xg1)
eta2_dp <- beta_dp %*% t(Xg2)
mu1_dp <- exp(eta1_dp)
mu2_dp <- exp(eta2_dp)
delta_dp <- mu2_dp - mu1_dp


## ===========================================================================
## A3. GLMM CURVES
## ===========================================================================

curve_draws_glmm <- function(fit, Xg1, Xg2, nsim = 5000L) {
  bhat <- as.numeric(fixef(fit)$cond)
  V <- as.matrix(vcov(fit)$cond)
  B <- rmvn_safe(nsim, bhat, V)
  
  eta1 <- B %*% t(Xg1)
  eta2 <- B %*% t(Xg2)
  
  mu1 <- exp(eta1)
  mu2 <- exp(eta2)
  
  list(mu1 = mu1, mu2 = mu2, delta = mu2 - mu1)
}

set.seed(1001)
cur_pois <- curve_draws_glmm(fit_pois, Xg1, Xg2, n_beta_sim)
set.seed(1002)
cur_nb <- curve_draws_glmm(fit_nb, Xg1, Xg2, n_beta_sim)

curve_table <- function(draws, model) {
  data.frame(
    model = model,
    day = day_grid,
    stage1_mean = colMeans(draws$mu1),
    stage1_lower = apply(draws$mu1, 2, quantile, .025),
    stage1_upper = apply(draws$mu1, 2, quantile, .975),
    stage2_mean = colMeans(draws$mu2),
    stage2_lower = apply(draws$mu2, 2, quantile, .025),
    stage2_upper = apply(draws$mu2, 2, quantile, .975),
    delta_mean = colMeans(draws$delta),
    delta_lower = apply(draws$delta, 2, quantile, .025),
    delta_upper = apply(draws$delta, 2, quantile, .975)
  )
}

curves_all <- bind_rows(
  curve_table(list(mu1 = mu1_dp, mu2 = mu2_dp, delta = delta_dp), "fcDPGLM"),
  curve_table(cur_pois, "Poisson GLMM"),
  curve_table(cur_nb, "NB GLMM")
)

write.csv(curves_all,
          file.path(out_dir, "full_data_stage_curves_and_delta.csv"),
          row.names = FALSE)


## ===========================================================================
## A4. STAGE CONTRASTS
## ===========================================================================

contrast_table <- function(draws, model) {
  avg_delta <- rowMeans(draws$delta)
  avg_log_rr <- rowMeans(log(draws$mu2) - log(draws$mu1))
  avg_rr <- exp(avg_log_rr)
  
  out <- rbind(
    average_delta = summ(avg_delta),
    average_log_rate_ratio = summ(avg_log_rr),
    average_rate_ratio = summ(avg_rr)
  )
  data.frame(model = model, contrast = rownames(out), out, row.names = NULL)
}

contrast_all <- bind_rows(
  contrast_table(list(mu1 = mu1_dp, mu2 = mu2_dp, delta = delta_dp), "fcDPGLM"),
  contrast_table(cur_pois, "Poisson GLMM"),
  contrast_table(cur_nb, "NB GLMM")
)

write.csv(contrast_all,
          file.path(out_dir, "full_data_scientific_contrasts.csv"),
          row.names = FALSE)

print(contrast_all)


## ===========================================================================
## A5. PLOTS: m_d(t) AND Delta(t)
## ===========================================================================

stage_long <- bind_rows(
  curves_all %>%
    transmute(model, day, stage = "Stage 1",
              mean = stage1_mean, lower = stage1_lower, upper = stage1_upper),
  curves_all %>%
    transmute(model, day, stage = "Stage 2",
              mean = stage2_mean, lower = stage2_lower, upper = stage2_upper)
)

p_curves <- ggplot(stage_long,
                   aes(day, mean, color = model, fill = model,
                       linetype = model, group = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .10,
              color = NA) +
  geom_line(linewidth = .9) +
  facet_wrap(~stage) +
  labs(
    x = "Day since start of stage",
    y = "Mean seizure count",
    color = "Model",
    fill = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave(file.path(out_dir, "compare_stage_curves.pdf"),
       p_curves, width = 8, height = 4.8)

p_delta <- ggplot(curves_all,
                  aes(day, delta_mean, color = model, fill = model,
                      linetype = model, group = model)) +
  geom_ribbon(aes(ymin = delta_lower, ymax = delta_upper),
              alpha = .10, color = NA) +
  geom_line(linewidth = .9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Day since start of stage",
    y = expression(Delta(t)),
    color = "Model",
    fill = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave(file.path(out_dir, "compare_delta.pdf"),
       p_delta, width = 7.5, height = 4.8)


## ===========================================================================
## A6. PREDICTIVE REPRODUCTION: ZEROS, VARIANCE, UPPER TAIL
## ===========================================================================

# DPGLM posterior predictive draws.
support <- fit_dp$support
theta_dp <- fit_dp$theta[keep, , drop = FALSE]
w_dp <- fit_dp$weights[keep, , drop = FALSE]

set.seed(2001)
pp_idx <- sample(seq_along(keep), min(n_pp_sim, length(keep)))
yrep_dp <- matrix(NA_integer_, length(pp_idx), nrow(dat))

for (j in seq_along(pp_idx)) {
  s <- pp_idx[j]
  prob <- tilted_pmf(theta_dp[s, ], w_dp[s, ], support)
  yrep_dp[j, ] <- draw_counts(prob, support)
}

# ---------------------------------------------------------------------------
# GLMM predictive draws given the fitted random effects.
#
# The checks concern the observed rats, so each GLMM's predicted rat effects
# are held fixed, as the fcDPGLM draws condition on the sampled rat effects.
# ---------------------------------------------------------------------------

mu_pois <- predict(
  fit_pois,
  type = "response",
  re.form = NULL
)

mu_nb <- predict(
  fit_nb,
  type = "response",
  re.form = NULL
)

# Poisson predictive draws.
set.seed(2002)
yrep_pois <- matrix(
  rpois(
    n_pp_sim * length(mu_pois),
    lambda = rep(mu_pois, each = n_pp_sim)
  ),
  nrow = n_pp_sim,
  ncol = length(mu_pois)
)

# NB2 predictive draws.
# For glmmTMB::nbinom2:
#   Var(Y_i | b_r) = mu_i + mu_i^2 / phi,
# where sigma(fit_nb) is phi.
phi_nb <- sigma(fit_nb)

set.seed(2003)
yrep_nb <- matrix(
  rnbinom(
    n_pp_sim * length(mu_nb),
    mu = rep(mu_nb, each = n_pp_sim),
    size = phi_nb
  ),
  nrow = n_pp_sim,
  ncol = length(mu_nb)
)

pp_all <- bind_rows(
  cbind(model = "fcDPGLM", pp_summary(dat$y, yrep_dp)),
  cbind(model = "Poisson GLMM", pp_summary(dat$y, yrep_pois)),
  cbind(model = "NB GLMM", pp_summary(dat$y, yrep_nb))
)

write.csv(pp_all,
          file.path(out_dir, "full_data_predictive_checks.csv"),
          row.names = FALSE)

cat("\nPredictive reproduction of distributional summaries:\n")
print(pp_all)


## ===========================================================================
## A7. PREDICTIVE REPRODUCTION RATIO PLOT
## ===========================================================================

# Put statistics with very different natural scales on a common scale:
#
#       predictive statistic / observed statistic.
#
# The mean is omitted; m_d(t) and Delta(t) show the mean comparison.

pp_ratio <- pp_all %>%
  mutate(
    ratio_median = predictive_median / observed,
    ratio_lower = lower95 / observed,
    ratio_upper = upper95 / observed,
    statistic = factor(
      statistic,
      levels = c("maximum", "p95", "p90", "variance", "zero_rate"),
      labels = c(
        "Maximum",
        "95th percentile",
        "90th percentile",
        "Variance",
        "Zero proportion"
    )
    )
  ) %>%
  filter(!is.na(statistic))

p_pp_ratio <- ggplot(
  pp_ratio,
  aes(
    x = ratio_median,
    y = statistic,
    color = model,
    group = model
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = ratio_lower,
      xmax = ratio_upper
    ),
    orientation = "y",
    width = 0.15,
    position = position_dodge(width = 0.55)
  ) +
  geom_point(
    size = 2.5,
    position = position_dodge(width = 0.55)
  ) +
  labs(
    x = "Predictive-to-observed ratio",
    y = NULL,
    color = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave(
  file.path(out_dir, "compare_predictive_reproduction.pdf"),
  p_pp_ratio,
  width = 7.5,
  height = 4.8
)


## ===========================================================================
## A8. STAGE-SPECIFIC PREDICTIVE CHECKS
## ===========================================================================

stage_pp <- function(yrep, model) {
  bind_rows(lapply(sort(unique(dat$stage)), function(s) {
    ii <- which(dat$stage == s)
    z <- pp_summary(dat$y[ii], yrep[, ii, drop = FALSE])
    cbind(model = model, stage = s, z)
  }))
}

pp_stage_all <- bind_rows(
  stage_pp(yrep_dp, "fcDPGLM"),
  stage_pp(yrep_pois, "Poisson GLMM"),
  stage_pp(yrep_nb, "NB GLMM")
)

write.csv(pp_stage_all,
          file.path(out_dir, "full_data_predictive_checks_by_stage.csv"),
          row.names = FALSE)



cat("\nDone: full-data comparison only. fcDPGLM was NOT refitted.\n")
cat("Results written to:", out_dir, "\n")
cat("  full_data_glmms.rds\n")
cat("  full_data_stage_curves_and_delta.csv\n")
cat("  full_data_scientific_contrasts.csv\n")
cat("  full_data_predictive_checks.csv\n")
cat("  full_data_predictive_checks_by_stage.csv\n")
cat("  compare_stage_curves.pdf\n")
cat("  compare_delta.pdf\n")
cat("  compare_predictive_reproduction.pdf\n")
