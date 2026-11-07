#!/usr/bin/env Rscript

# Posterior summaries for the fcDPGLM fit: regression coefficients, stage
# curves m_d(t), the stage difference Delta(t), rate ratios, the shape
# contrast S(t), rat effects, posterior predictive checks and MCMC traces.
# Every quantity is a functional of the posterior draws in fit_path.
#
# Usage (from the repository root, after scripts/01_fit_fcdpglm.R):
#   Rscript scripts/03_posterior_summaries.R [fit.rds] [out_dir]
#
# Model:
# E(Y_i | d_i, t_i, e_i, a_r)
#   = e_i * exp(
#       b0
#       + b_stage * 1{d=2}
#       + B(t)'g1 * 1{d=1}
#       + B(t)'g2 * 1{d=2}
#       + a_r
#     )
#
# Because theta is solved so that b'(theta_i) = mu_i, every mean
# functional depends only on (beta, a); the CRM weights matter only
# for the predictive distribution of counts.


## ========================================================================
## SETTINGS
## ========================================================================

fit_path   <- "results/fcdpglm_fit.rds"
out_dir    <- "results/fcdpglm"
# Optional command-line override: <fit.rds> <out_dir>
cl_args <- commandArgs(trailingOnly = TRUE)
if (length(cl_args) >= 1L) fit_path <- cl_args[1]
if (length(cl_args) >= 2L) out_dir  <- cl_args[2]
burn_in    <- NULL     # NULL = use fit$tuning$adapt_until; integer overrides
n_pp_draws <- 200L     # posterior-predictive replicate datasets

# Save all fcDPGLM figures and result tables here (out_dir set above).
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260921)

# Tables print to the console; figures are saved to out_dir.  When run with
# Rscript, on-screen plots are discarded.
if (!interactive()) grDevices::pdf(NULL)

suppressPackageStartupMessages({
  library(splines)
  library(ggplot2)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)


## ========================================================================
## LOAD FIT
## ========================================================================

fit_object <- readRDS(fit_path)

fit <- fit_object$rat_dpglm
dat <- as.data.frame(fit_object$data)
X   <- as.matrix(fit_object$design_matrix)
cn  <- colnames(X)

n       <- nrow(dat)
support <- fit$support


## ========================================================================
## DESIGN MATRIX INDICES
## ========================================================================

idx_int   <- which(cn == "intercept")
idx_stage <- which(cn == "stage2")
idx_g1    <- grep("^stage1_spline", cn)
idx_g2    <- grep("^stage2_spline", cn)

stopifnot(
  length(idx_int) == 1L,
  length(idx_stage) == 1L,
  length(idx_g1) == length(idx_g2),
  length(idx_g1) > 0L
)


## ========================================================================
## POSTERIOR DRAWS
## ========================================================================

n_iter <- nrow(fit$beta)

# Adaptation ends at fit$tuning$adapt_until (0 for the pooled fit, whose
# draws are all post-adaptation).
if (is.null(burn_in)) {
  if (!is.null(fit$tuning$adapt_until)) {
    burn_in <- as.integer(fit$tuning$adapt_until)
  } else {
    burn_in <- floor(n_iter / 2)
  }
}

if (burn_in < 0L || burn_in >= n_iter)
  stop("burn_in must be between 0 and n_iter - 1.", call. = FALSE)

keep <- seq.int(burn_in + 1L, n_iter)
S    <- length(keep)

beta_draws  <- fit$beta[keep, , drop = FALSE]
a_draws     <- fit$a[keep, , drop = FALSE]
theta_draws <- fit$theta[keep, , drop = FALSE]
w_draws     <- fit$weights[keep, , drop = FALSE]

colnames(beta_draws) <- cn

# fit$a columns follow unique(rat_id) = order of first appearance
rat_levels <- unique(as.character(dat$int_id))
colnames(a_draws) <- rat_levels


## ========================================================================
## UTILITIES
## ========================================================================

summ <- function(x) {
  
  q <- unname(
    stats::quantile(
      x,
      c(0.025, 0.5, 0.975)
    )
  )
  
  c(
    mean    = mean(x),
    median  = q[2],
    lower95 = q[1],
    upper95 = q[3],
    sd      = stats::sd(x)
  )
}


summ_cols <- function(M) {
  as.data.frame(
    t(apply(M, 2L, summ))
  )
}


# Geyer initial-positive-sequence effective sample size
ess_basic <- function(x) {
  x  <- as.numeric(x)
  nx <- length(x)
  
  if (nx < 10L || stats::sd(x) == 0)
    return(NA_real_)
  
  a <- stats::acf(
    x,
    lag.max = min(500L, nx - 2L),
    plot = FALSE
  )$acf[-1L]
  
  total <- 0
  
  for (k in seq(1L, length(a) - 1L, by = 2L)) {
    
    pair <- a[k] + a[k + 1L]
    
    if (pair < 0)
      break
    
    total <- total + pair
  }
  nx / (1 + 2 * total)
}


# Simultaneous max-|z| band for draws x grid matrix
sim_band <- function(M, level = 0.95) {
  
  m <- colMeans(M)
  s <- apply(M, 2L, stats::sd)
  
  s[s <= 0] <- .Machine$double.eps
  
  z <- apply(
    abs(sweep(M, 2L, m, "-")) /
      rep(s, each = nrow(M)),
    1L,
    max
  )
  
  q <- unname(
    stats::quantile(z, level)
  )
  
  list(
    mean  = m,
    lower = m - q * s,
    upper = m + q * s
  )
}


# One count per row from matrix of row-wise probabilities
draw_counts <- function(prob, values) {
  
  cs <- t(
    apply(prob, 1L, cumsum)
  )
  
  cs[, ncol(cs)] <- 1
  
  values[
    max.col(
      cs >= stats::runif(nrow(prob)),
      ties.method = "first"
    )
  ]
}


# Row-wise tilted pmf:
# p(k) proportional to w_k exp(theta*k)
tilted_pmf <- function(theta, weights) {
  
  lp <- outer(theta, support) +
    rep(log(weights), each = length(theta))
  
  lp <- lp - apply(lp, 1L, max)
  
  p <- exp(lp)
  
  p / rowSums(p)
}


## ========================================================================
## 0. SAMPLER DIAGNOSTICS
## ========================================================================

cat(
  "\n================ 0. SAMPLER DIAGNOSTICS ================\n"
)

cat(
  sprintf(
    "iterations = %d, burn-in = %d, draws used = %d\n",
    n_iter,
    burn_in,
    S
  )
)

if (!is.null(fit$acceptance$post_adaptation$joint) &&
    is.finite(fit$acceptance$post_adaptation$joint)) {

  cat(
    sprintf(
      paste0(
        "post-adaptation acceptance: joint(beta,a) = %.3f | ",
        "CRM = %.3f\n"
      ),
      fit$acceptance$post_adaptation$joint,
      fit$acceptance$post_adaptation$crm
    )
  )

} else {

  cat(
    sprintf(
      paste0(
        "acceptance: beta = %.3f | ",
        "random effects = %.3f | CRM = %.3f\n"
      ),
      fit$acceptance$beta,
      fit$acceptance$random_effects,
      fit$acceptance$crm
    )
  )
}


crm_mass <- rowSums(fit$weights)

cat(
  sprintf(
    paste0(
      "CRM total mass: start = %.1f, ",
      "end of adaptation = %.2f, final = %.2f\n"
    ),
    crm_mass[1],
    crm_mass[max(1L, burn_in)],
    crm_mass[n_iter]
  )
)

cat(
  paste0(
    "CRM total mass is a scale direction of the normalized likelihood; ",
    "use it only as a gross stability diagnostic.\n",
    "Convergence assessment below focuses on beta and scientifically ",
    "relevant posterior functionals.\n"
  )
)

# Scale-free CRM weights pi_k = w_k / sum_j w_j (invariant to rescaling w).
w_norm <- sweep(
  w_draws,
  1L,
  rowSums(w_draws),
  "/"
)

crm_norm_mean <- colMeans(w_norm)
crm_norm_ess  <- apply(w_norm, 2L, ess_basic)

# Report only support points carrying non-negligible posterior mean mass.
crm_active <- which(crm_norm_mean >= 0.01)

if (length(crm_active) > 0L) {
  crm_norm_tab <- data.frame(
    support = support[crm_active],
    posterior_mean_weight = crm_norm_mean[crm_active],
    ESS = crm_norm_ess[crm_active]
  )

  cat("\nNormalized CRM weights with posterior mean >= 0.01:\n")
  print(
    transform(
      crm_norm_tab,
      posterior_mean_weight = round(posterior_mean_weight, 4),
      ESS = round(ESS, 1)
    ),
    row.names = FALSE
  )
}


ess_beta <- apply(
  beta_draws,
  2L,
  ess_basic
)

print(
  round(
    data.frame(
      ESS = ess_beta,
      ESS_per_1000 = 1000 * ess_beta / S
    ),
    1
  )
)

# ESS and MCSE of the mean functionals are computed further below.


## ------------------------------------------------------------------------
## Diagnostic plots
## ------------------------------------------------------------------------

diag_df <- data.frame(
  iteration = seq_len(n_iter),
  loglik     = fit$loglik,
  crm_mass   = crm_mass,
  crm_log_r  = fit$crm_log_r
)


# Log likelihood

p_loglik <- ggplot(
  diag_df,
  aes(x = iteration, y = loglik)
) +
  geom_line(linewidth = 0.4) +
  geom_vline(
    xintercept = burn_in,
    linetype = "dashed"
  ) +
  labs(
    x = "Iteration",
    y = "Log likelihood",
    title = "Log likelihood"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_loglik)


# CRM total mass

p_crm_mass <- ggplot(
  diag_df,
  aes(x = iteration, y = crm_mass)
) +
  geom_line(linewidth = 0.4) +
  geom_vline(
    xintercept = burn_in,
    linetype = "dashed"
  ) +
  scale_y_log10() +
  labs(
    x = "Iteration",
    y = expression(sum(w)),
    title = "CRM total mass",
    subtitle = "Scale direction; not a primary convergence diagnostic"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_crm_mass)


# Scale-free CRM composition: trace the support points with the largest
# posterior mean normalized weights.  One panel per support point.
n_crm_trace <- min(5L, ncol(w_norm))
crm_top <- order(crm_norm_mean, decreasing = TRUE)[seq_len(n_crm_trace)]

crm_norm_plots <- lapply(
  crm_top,
  function(k) {

    df_k <- data.frame(
      draw = seq_len(S),
      normalized_weight = w_norm[, k]
    )

    ggplot(
      df_k,
      aes(x = draw, y = normalized_weight)
    ) +
      geom_line(linewidth = 0.35) +
      labs(
        x = "Retained draw",
        y = expression(pi[k]),
        title = bquote(
          pi[.(support[k])] ==
            w[.(support[k])] / sum(j, w[j])
        ),
        subtitle = sprintf(
          "Posterior mean = %.3f; ESS = %.1f",
          crm_norm_mean[k],
          crm_norm_ess[k]
        )
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank()
      )
  }
)

# One plot per support point.
for (p in crm_norm_plots)
  print(p)

# All traces in one grid.
if (has_patchwork) {

  p_crm_norm_grid <- patchwork::wrap_plots(
    crm_norm_plots,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = "Leading normalized CRM-weight traces"
    )

  print(p_crm_norm_grid)

} else {

  cat(
    paste0(
      "\\nPackage 'patchwork' is not installed, so the five normalized ",
      "CRM traces were printed separately only.\\n",
      "Install patchwork if you also want the combined diagnostic grid.\\n"
    )
  )
}




# CRM log acceptance ratio

p_crm_accept <- ggplot(
  diag_df,
  aes(x = iteration, y = crm_log_r)
) +
  geom_line(linewidth = 0.4) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_cartesian(
    ylim = c(-100, 20)
  ) +
  labs(
    x = "Iteration",
    y = "Log R",
    title = "CRM log acceptance ratio"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_crm_accept)


# Trace plots: intercept and stage effect

trace_df <- rbind(
  
  data.frame(
    draw = seq_len(S),
    value = beta_draws[, idx_int],
    parameter = "Intercept"
  ),
  
  data.frame(
    draw = seq_len(S),
    value = beta_draws[, idx_stage],
    parameter = "Stage 2"
  )
)


p_trace <- ggplot(
  trace_df,
  aes(
    x = draw,
    y = value,
    group = parameter,
    linetype = parameter
  )
) +
  geom_line(linewidth = 0.45) +
  labs(
    x = "Retained draw",
    y = "Coefficient",
    linetype = NULL,
    title = "Trace plots"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(p_trace)


## ========================================================================
## 1. COEFFICIENT POSTERIOR SUMMARIES
## ========================================================================

cat(
  "\n================ 1. COEFFICIENTS ================\n"
)

beta_tab <- summ_cols(beta_draws)

beta_tab$`P(>0)` <- colMeans(
  beta_draws > 0
)

print(
  round(beta_tab, 4)
)

write.csv(
  cbind(parameter = rownames(beta_tab), beta_tab),
  file.path(out_dir, "coefficient_posterior_summaries.csv"),
  row.names = FALSE
)


# Priors are those of scripts/01_fit_fcdpglm.R; the preliminary NB GLMM
# provides starting values and the G0 size, not a prior.


## ========================================================================
## 2. FITTED STAGE-SPECIFIC CURVES f1(t), f2(t)
## ========================================================================

# t was standardized within stage, so compare stages at a common
# calendar day (days since start of that stage) and map each day
# to its stage's t.

basis <- ns(
  dat$t,
  df = 3,
  intercept = FALSE
)


X_rebuilt <- cbind(
  1,
  dat$stage2,
  predict(basis, dat$t) * (dat$stage == 1L),
  predict(basis, dat$t) * (dat$stage == 2L)
)

stopifnot(
  max(abs(X_rebuilt - X)) < 1e-8
)


t_of_day <- function(day, s) {
  
  d <- dat$day[
    dat$stage == s
  ]
  
  (day - mean(d)) / stats::sd(d)
}


rng <- sapply(
  split(dat$day, dat$stage),
  range
)

day_grid <- seq(
  max(rng[1, ]),
  min(rng[2, ])
)

G <- length(day_grid)


B1 <- predict(
  basis,
  t_of_day(day_grid, 1L)
)

B2 <- predict(
  basis,
  t_of_day(day_grid, 2L)
)


# Linear predictors for one full day exposure and typical rat (a = 0)

eta1 <-
  beta_draws[, idx_g1, drop = FALSE] %*% t(B1) +
  beta_draws[, idx_int]

eta2 <-
  beta_draws[, idx_g2, drop = FALSE] %*% t(B2) +
  beta_draws[, idx_int] +
  beta_draws[, idx_stage]


mu1 <- exp(eta1)
mu2 <- exp(eta2)


# Population-average versions marginal over a ~ N(0, sigma_a^2).
# If sigma_a was sampled, the factor exp(sigma_a^2 / 2) is applied draw by draw.

sigma_a_sampled <- !is.null(fit$sigma_a_samples) &&
  length(unique(fit$sigma_a_samples[keep])) > 1L
sa_draws  <- if (sigma_a_sampled) fit$sigma_a_samples[keep] else rep(fit$sigma_a, S)
infl_draw <- exp(sa_draws^2 / 2)
infl      <- mean(infl_draw)


curve_tab <- data.frame(
  
  day = day_grid,
  
  stage1_mean =
    colMeans(mu1),
  
  stage1_lower =
    apply(mu1, 2L, stats::quantile, 0.025),
  
  stage1_upper =
    apply(mu1, 2L, stats::quantile, 0.975),
  
  stage2_mean =
    colMeans(mu2),
  
  stage2_lower =
    apply(mu2, 2L, stats::quantile, 0.025),
  
  stage2_upper =
    apply(mu2, 2L, stats::quantile, 0.975),
  
  stage1_marginal =
    colMeans(mu1 * infl_draw),
  
  stage2_marginal =
    colMeans(mu2 * infl_draw)
)


cat(
  "\n================ 2. STAGE CURVES ",
  "(counts/day, typical rat) ================\n"
)


print(
  round(
    curve_tab[
      round(seq(1, G, length.out = 8)),
      1:7
    ],
    3
  ),
  row.names = FALSE
)

write.csv(
  curve_tab,
  file.path(out_dir, "stage_specific_curves.csv"),
  row.names = FALSE
)


## ------------------------------------------------------------------------
## Stage-specific curve plot
## ------------------------------------------------------------------------

curve_plot_df <- rbind(
  
  data.frame(
    day   = day_grid,
    mean  = colMeans(mu1),
    lower = apply(
      mu1,
      2L,
      stats::quantile,
      0.025
    ),
    upper = apply(
      mu1,
      2L,
      stats::quantile,
      0.975
    ),
    stage = "Stage 1"
  ),
  
  data.frame(
    day   = day_grid,
    mean  = colMeans(mu2),
    lower = apply(
      mu2,
      2L,
      stats::quantile,
      0.025
    ),
    upper = apply(
      mu2,
      2L,
      stats::quantile,
      0.975
    ),
    stage = "Stage 2"
  )
)


p_stage <- ggplot(
  curve_plot_df,
  aes(
    x = day,
    y = mean,
    group = stage,
    color = stage,
    fill = stage
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  labs(
    x = "Day since start of stage",
    y = "Mean seizure count",
    color = NULL,
    fill = NULL
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(p_stage)
ggsave(file.path(out_dir, "stage_specific_curves.pdf"),
       p_stage, width = 7.5, height = 4.8)


## ========================================================================
## 3. DELTA(t)
## ========================================================================

cat(
  "\n================ 3. DELTA(t) ================\n"
)


# Rate difference: counts/day
delta_mu <- mu2 - mu1

# Log rate ratio
log_rr <- eta2 - eta1


delta_tab <- data.frame(
  
  day = day_grid,
  
  delta_mean =
    colMeans(delta_mu),
  
  delta_median =
    apply(
      delta_mu,
      2L,
      stats::median
    ),
  
  delta_lower =
    apply(
      delta_mu,
      2L,
      stats::quantile,
      0.025
    ),
  
  delta_upper =
    apply(
      delta_mu,
      2L,
      stats::quantile,
      0.975
    ),
  
  `P(delta<0)` =
    colMeans(delta_mu < 0),
  
  rate_ratio =
    exp(colMeans(log_rr)),
  
  rr_lower =
    exp(
      apply(
        log_rr,
        2L,
        stats::quantile,
        0.025
      )
    ),
  
  rr_upper =
    exp(
      apply(
        log_rr,
        2L,
        stats::quantile,
        0.975
      )
    ),
  
  check.names = FALSE
)


print(
  round(
    delta_tab[
      round(seq(1, G, length.out = 8)),
    ],
    3
  ),
  row.names = FALSE
)

write.csv(
  delta_tab,
  file.path(out_dir, "delta_by_day.csv"),
  row.names = FALSE
)


cat(
  "\nMarginal (population-average) Delta = conditional Delta x",
  round(infl, 3),
  "\n"
)


## ------------------------------------------------------------------------
## Delta(t) plot
## ------------------------------------------------------------------------

sb <- sim_band(
  delta_mu
)


delta_plot_df <- data.frame(
  
  day = day_grid,
  
  mean = sb$mean,
  
  point_lower =
    apply(
      delta_mu,
      2L,
      stats::quantile,
      0.025
    ),
  
  point_upper =
    apply(
      delta_mu,
      2L,
      stats::quantile,
      0.975
    ),
  
  simultaneous_lower =
    sb$lower,
  
  simultaneous_upper =
    sb$upper
)


p_delta <- ggplot(
  delta_plot_df,
  aes(
    x = day,
    y = mean
  )
) +
  geom_ribbon(
    aes(
      ymin = simultaneous_lower,
      ymax = simultaneous_upper,
      fill = "Simultaneous 95%"
    ),
    alpha = 0.15
  ) +
  geom_ribbon(
    aes(
      ymin = point_lower,
      ymax = point_upper,
      fill = "Pointwise 95%"
    ),
    alpha = 0.30
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Day since start of stage",
    y = expression(Delta(t)),
    fill = NULL
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(p_delta)
ggsave(file.path(out_dir, "delta_curve.pdf"),
       p_delta, width = 7.5, height = 4.8)


## ========================================================================
## 4. AVERAGE STAGE DIFFERENCE AND DIRECTION
## ========================================================================

cat(
  "\n========== 4. AVERAGE STAGE DIFFERENCE / ",
  "4. DIRECTION ==========\n"
)


# (a) coefficient itself

b_stage <- beta_draws[, idx_stage]


# (b) average over common day grid

avg_log_rr <- rowMeans(
  log_rr
)

avg_delta <- rowMeans(
  delta_mu
)


# (c) Standardised over observed design:
# every row evaluated under both stages at own t,
# typical rat, one day exposure.

Xall <- predict(
  basis,
  dat$t
)


eta_s1_all <-
  beta_draws[, idx_g1, drop = FALSE] %*% t(Xall) +
  beta_draws[, idx_int]


eta_s2_all <-
  beta_draws[, idx_g2, drop = FALSE] %*% t(Xall) +
  beta_draws[, idx_int] +
  beta_draws[, idx_stage]


gcomp_delta <-
  rowMeans(exp(eta_s2_all)) -
  rowMeans(exp(eta_s1_all))


## ------------------------------------------------------------------------
## Functional convergence diagnostics
## ------------------------------------------------------------------------
# Computed from the retained draws.

stopifnot(
  length(avg_delta) == S,
  length(gcomp_delta) == S,
  length(avg_log_rr) == S,
  nrow(delta_mu) == S,
  ncol(delta_mu) == length(day_grid)
)

functional_ess <- c(
  avg_delta       = ess_basic(avg_delta),
  gcomp_delta     = ess_basic(gcomp_delta),
  avg_log_rr      = ess_basic(avg_log_rr),
  avg_rate_ratio  = ess_basic(exp(avg_log_rr))
)

delta_ess <- apply(
  delta_mu,
  2L,
  ess_basic
)

functional_sd <- c(
  avg_delta       = stats::sd(avg_delta),
  gcomp_delta     = stats::sd(gcomp_delta),
  avg_log_rr      = stats::sd(avg_log_rr),
  avg_rate_ratio  = stats::sd(exp(avg_log_rr))
)

functional_mcse <- functional_sd / sqrt(functional_ess)

functional_diag <- data.frame(
  ESS = functional_ess,
  ESS_per_1000 = 1000 * functional_ess / S,
  posterior_SD = functional_sd,
  MCSE_mean = functional_mcse
)

cat("\nESS / MCSE for primary inferential functionals:\n")
print(round(functional_diag, 3))

cat("\nESS for Delta(t):\n")
print(
  data.frame(
    day = day_grid,
    ESS = round(delta_ess, 1)
  ),
  row.names = FALSE
)


# Traces of the main functionals
functional_trace_df <- rbind(
  data.frame(
    draw = seq_len(S),
    value = avg_delta,
    functional = "Average Delta"
  ),
  data.frame(
    draw = seq_len(S),
    value = gcomp_delta,
    functional = "G-computation Delta"
  ),
  data.frame(
    draw = seq_len(S),
    value = avg_log_rr,
    functional = "Average log rate ratio"
  )
)

p_functional_trace <- ggplot(
  functional_trace_df,
  aes(x = draw, y = value)
) +
  geom_line(linewidth = 0.30) +
  facet_wrap(
    ~ functional,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    x = "Retained draw",
    y = NULL,
    title = "MCMC traces for primary inferential functionals"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

print(p_functional_trace)


# Pointwise Delta(t) traces at early, middle, and late follow-up.
delta_trace_days <- unique(
  c(
    min(day_grid),
    day_grid[which.min(abs(day_grid - stats::median(day_grid)))],
    max(day_grid)
  )
)

delta_trace_idx <- match(delta_trace_days, day_grid)

delta_trace_df <- do.call(
  rbind,
  lapply(
    seq_along(delta_trace_idx),
    function(j) {
      data.frame(
        draw = seq_len(S),
        value = delta_mu[, delta_trace_idx[j]],
        day = factor(
          paste0("Day ", delta_trace_days[j]),
          levels = paste0("Day ", delta_trace_days)
        )
      )
    }
  )
)

p_delta_trace <- ggplot(
  delta_trace_df,
  aes(x = draw, y = value)
) +
  geom_line(linewidth = 0.30) +
  facet_wrap(
    ~ day,
    scales = "free_y",
    ncol = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  labs(
    x = "Retained draw",
    y = expression(Delta(t)),
    title = expression("MCMC traces for pointwise " * Delta(t))
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

print(p_delta_trace)


# Autocorrelation of two mean functionals.
acf_to_df <- function(x, label, lag_max = 100L) {
  ac <- stats::acf(
    x,
    lag.max = lag_max,
    plot = FALSE,
    demean = TRUE
  )

  data.frame(
    lag = as.numeric(ac$lag),
    acf = as.numeric(ac$acf),
    functional = label
  )
}

acf_df <- rbind(
  acf_to_df(avg_delta, "Average Delta"),
  acf_to_df(avg_log_rr, "Average log rate ratio")
)

p_functional_acf <- ggplot(
  acf_df,
  aes(x = lag, y = acf)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  geom_segment(
    aes(
      xend = lag,
      y = 0,
      yend = acf
    ),
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ functional,
    ncol = 1
  ) +
  labs(
    x = "Lag",
    y = "Autocorrelation",
    title = "Autocorrelation of primary inferential functionals"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

print(p_functional_acf)


q1 <- rbind(
  
  `beta_stage (log scale)` =
    summ(b_stage),
  
  `mean log rate ratio over days` =
    summ(avg_log_rr),
  
  `mean rate ratio over days` =
    summ(exp(avg_log_rr)),
  
  `mean Delta over days (cts/day)` =
    summ(avg_delta),
  
  `g-computation Delta (cts/day)` =
    summ(gcomp_delta)
)


q1 <- cbind(
  
  q1,
  
  `P(<0)` = c(
    
    mean(b_stage < 0),
    
    mean(avg_log_rr < 0),
    
    mean(exp(avg_log_rr) < 1),
    
    mean(avg_delta < 0),
    
    mean(gcomp_delta < 0)
  )
)


print(
  round(q1, 4)
)

write.csv(
  cbind(contrast = rownames(q1), q1),
  file.path(out_dir, "scientific_contrasts.csv"),
  row.names = FALSE
)


prob_report <- function(indicator, label) {

  n_event <- sum(indicator)
  p_hat <- mean(indicator)

  if (n_event == 0L) {
    cat(
      sprintf(
        "  %-42s no events in %d retained draws (Monte Carlo estimate = 0)\n",
        label,
        length(indicator)
      )
    )
  } else {
    cat(
      sprintf(
        "  %-42s %.4f (%d / %d draws)\n",
        label,
        p_hat,
        n_event,
        length(indicator)
      )
    )
  }
}

cat(
  "\nPosterior draw frequencies for Stage 2 being LOWER than Stage 1:\n"
)

prob_report(
  b_stage < 0,
  "beta_stage < 0:"
)

prob_report(
  avg_delta < 0,
  "average Delta < 0:"
)

prob_report(
  apply(delta_mu < 0, 1L, all),
  "Delta(t) < 0 for every day:"
)

prob_report(
  apply(delta_mu < 0, 1L, any),
  "Delta(t) < 0 for at least one day:"
)



## ========================================================================
## 5. DO TIME PROFILES DIFFER ACROSS STAGES?
## ========================================================================

cat(
  "\n================ 5. TIME PROFILES ================\n"
)


d_gamma <-
  beta_draws[, idx_g2, drop = FALSE] -
  beta_draws[, idx_g1, drop = FALSE]


colnames(d_gamma) <- paste0(
  "gamma2_minus_gamma1_",
  seq_len(ncol(d_gamma))
)


print(
  round(
    cbind(
      summ_cols(d_gamma),
      `P(>0)` = colMeans(d_gamma > 0)
    ),
    4
  )
)


# Shape comparison: simultaneous band for the centered difference of the stage curves.


## ------------------------------------------------------------------------
## Shape-only comparison
## ------------------------------------------------------------------------

# Center each curve so only the time pattern is compared

shape1 <- eta1 - rowMeans(eta1)
shape2 <- eta2 - rowMeans(eta2)

shape_diff <- shape2 - shape1

sbs <- sim_band(
  shape_diff
)


cat(
  sprintf(
    paste0(
      "Simultaneous 95%% band for the shape difference ",
      "excludes zero somewhere: %s\n"
    ),
    if (
      any(
        sbs$lower > 0 |
        sbs$upper < 0
      )
    ) {
      "YES"
    } else {
      "no"
    }
  )
)


shape_plot_df <- data.frame(
  
  day = day_grid,
  
  mean = sbs$mean,
  
  lower = sbs$lower,
  
  upper = sbs$upper
)


p_shape <- ggplot(
  shape_plot_df,
  aes(
    x = day,
    y = mean
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Day since start of stage",
    y = "Centered log-rate difference"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_shape)
ggsave(file.path(out_dir, "time_profile_shape_difference.pdf"),
       p_shape, width = 7.5, height = 4.8)


## ========================================================================
## 6. BETWEEN-RAT HETEROGENEITY
## ========================================================================

cat(
  "\n================ 6. RAT HETEROGENEITY ================\n"
)


a_tab <- summ_cols(
  a_draws
)

a_tab$multiplier <- exp(
  a_tab$mean
)

print(
  round(a_tab, 3)
)


sd_a_draw <- apply(
  a_draws,
  1L,
  stats::sd
)


fold_range <- exp(
  apply(
    a_draws,
    1L,
    function(z) {
      max(z) - min(z)
    }
  )
)


het <- rbind(
  
  `realised SD of a_r` =
    summ(sd_a_draw),
  
  `fold range exp(max a - min a)` =
    summ(fold_range)
)


print(
  round(het, 3)
)

write.csv(
  cbind(summary = rownames(het), het),
  file.path(out_dir, "rat_heterogeneity_summary.csv"),
  row.names = FALSE
)


if (sigma_a_sampled) {
  q_sa <- stats::quantile(sa_draws, c(0.025, 0.5, 0.975))
  cat(sprintf(
    "\nsigma_a is SAMPLED: posterior median %.3f (95%% CI %.3f, %.3f); realised SD of the rat effects (posterior mean) %.3f.\n",
    q_sa[2], q_sa[1], q_sa[3], mean(sd_a_draw)))
} else {
cat(
  sprintf(
    paste0(
      "\nsigma_a is FIXED at %.4f ",
      "(empirical Bayes from the GLMM); ",
      "it is not sampled,\n"
    ),
    fit$sigma_a
  )
)


cat(
  "so the numbers above do not propagate uncertainty ",
  "in sigma_a itself.\n"
)


cat(
  sprintf(
    "Realised SD (posterior mean %.3f) vs fixed sigma_a %.3f: %s\n",
    
    mean(sd_a_draw),
    
    fit$sigma_a,
    
    if (
      mean(sd_a_draw) >
      1.3 * fit$sigma_a
    ) {
      
      paste0(
        "the data want MORE spread than the prior allows ",
        "- refit with a larger sigma_a"
      )
      
    } else if (
      mean(sd_a_draw) <
      0.7 * fit$sigma_a
    ) {
      
      "the prior is looser than the data need"
      
    } else {
      
      "consistent"
    }
  )
)
}


# Share of total variation in linear predictor attributable to rats

lp_fixed <- beta_draws %*% t(X)


var_share <-
  apply(
    a_draws,
    1L,
    stats::var
  ) /
  (
    apply(
      a_draws,
      1L,
      stats::var
    ) +
      apply(
        lp_fixed,
        1L,
        stats::var
      )
  )


cat(
  sprintf(
    paste0(
      "Share of linear-predictor variance from rat effects: ",
      "%.3f (95%% CI %.3f-%.3f)\n"
    ),
    
    mean(var_share),
    
    stats::quantile(
      var_share,
      0.025
    ),
    
    stats::quantile(
      var_share,
      0.975
    )
  )
)


## ------------------------------------------------------------------------
## Rat-specific effects: caterpillar plots
## ------------------------------------------------------------------------

# Posterior summaries on the model's log-rate scale.
rat_effect_summary <- data.frame(
  rat = rat_levels,
  mean = apply(a_draws, 2L, mean),
  median = apply(a_draws, 2L, stats::median),
  lower95 = apply(a_draws, 2L, stats::quantile, probs = 0.025),
  upper95 = apply(a_draws, 2L, stats::quantile, probs = 0.975),
  stringsAsFactors = FALSE
)

# Order rats from lowest to highest posterior mean.
rat_effect_summary <- rat_effect_summary[
  order(rat_effect_summary$mean),
  ,
  drop = FALSE
]

rat_effect_summary$rat_ordered <- factor(
  rat_effect_summary$rat,
  levels = rat_effect_summary$rat
)

p_rat_effects <- ggplot(
  rat_effect_summary,
  aes(x = mean, y = rat_ordered)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(xmin = lower95, xmax = upper95),
    width = 0,
    linewidth = 0.7,
    orientation = "y"
  ) +
  geom_point(size = 2.2) +
  labs(
    x = expression("Rat-specific random effect  " * a[r]),
    y = "Rat"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(p_rat_effects)
ggsave(file.path(out_dir, "rat_random_effects_log_scale.pdf"),
       p_rat_effects, width = 6, height = 6)


# Posterior summaries on the multiplicative rate scale exp(a_r).
# Posterior mean of exp(a_r), not exp of the posterior mean of a_r.
a_multiplier_draws <- exp(a_draws)

rat_multiplier_summary <- data.frame(
  rat = rat_levels,
  mean = apply(a_multiplier_draws, 2L, mean),
  median = apply(a_multiplier_draws, 2L, stats::median),
  lower95 = apply(
    a_multiplier_draws,
    2L,
    stats::quantile,
    probs = 0.025
  ),
  upper95 = apply(
    a_multiplier_draws,
    2L,
    stats::quantile,
    probs = 0.975
  ),
  stringsAsFactors = FALSE
)

rat_multiplier_summary <- rat_multiplier_summary[
  order(rat_multiplier_summary$mean),
  ,
  drop = FALSE
]

rat_multiplier_summary$rat_ordered <- factor(
  rat_multiplier_summary$rat,
  levels = rat_multiplier_summary$rat
)

cat("\nPosterior rat-specific multiplicative effects exp(a_r):\n")
rat_multiplier_print <- rat_multiplier_summary[
  ,
  c("rat", "mean", "median", "lower95", "upper95")
]
rat_multiplier_print[, -1] <- lapply(
  rat_multiplier_print[, -1, drop = FALSE],
  round,
  digits = 3
)
print(rat_multiplier_print, row.names = FALSE)


# exp(a_r) is the multiplicative effect of rat r on the expected rate.
p_rat_multiplier <- ggplot(
  rat_multiplier_summary,
  aes(x = mean, y = rat_ordered)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(xmin = lower95, xmax = upper95),
    width = 0,
    linewidth = 0.7,
    orientation = "y"
  ) +
  geom_point(size = 2.2) +
  labs(
    x = expression("Multiplicative rat effect  " * exp(a[r])),
    y = "Rat"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(p_rat_multiplier)
ggsave(file.path(out_dir, "rat_random_effects.pdf"),
       p_rat_multiplier, width = 6, height = 6)



## ========================================================================
## 7. POSTERIOR PREDICTIVE COUNTS
## ========================================================================

cat(
  "\n================ 7. POSTERIOR PREDICTIVE ================\n"
)

pp_idx <- sort(
  sample(
    seq_len(S),
    min(n_pp_draws, S)
  )
)

stat_names <- c(
  "mean",
  "variance",
  "zero_rate",
  "p90",
  "maximum"
)

pp <- matrix(
  NA_real_,
  length(pp_idx),
  length(stat_names),
  dimnames = list(
    NULL,
    stat_names
  )
)

yrep_store <- matrix(
  NA_real_,
  length(pp_idx),
  n
)

for (j in seq_along(pp_idx)) {
  
  s <- pp_idx[j]
  
  p <- tilted_pmf(
    theta_draws[s, ],
    w_draws[s, ]
  )
  
  yrep <- draw_counts(
    p,
    support
  )
  
  yrep_store[j, ] <- yrep
  
  pp[j, ] <- c(
    mean(yrep),
    stats::var(yrep),
    mean(yrep == 0),
    unname(
      stats::quantile(
        yrep,
        0.9
      )
    ),
    max(yrep)
  )
}


## ------------------------------------------------------------------------
## Posterior predictive summaries
## ------------------------------------------------------------------------

obs <- c(
  mean(dat$y),
  stats::var(dat$y),
  mean(dat$y == 0),
  unname(
    stats::quantile(
      dat$y,
      0.9
    )
  ),
  max(dat$y)
)

pp_tab <- data.frame(
  statistic = stat_names,
  
  observed = obs,
  
  post_median =
    apply(
      pp,
      2L,
      stats::median
    ),
  
  lower95 =
    apply(
      pp,
      2L,
      stats::quantile,
      0.025
    ),
  
  upper95 =
    apply(
      pp,
      2L,
      stats::quantile,
      0.975
    ),
  
  bayes_p =
    colMeans(
      sweep(
        pp,
        2L,
        obs,
        ">="
      )
    ),
  
  check.names = FALSE
)


# Round numeric columns only
pp_tab_print <- pp_tab

num_cols <- vapply(
  pp_tab_print,
  is.numeric,
  logical(1)
)

pp_tab_print[num_cols] <- lapply(
  pp_tab_print[num_cols],
  round,
  digits = 3
)

print(
  pp_tab_print,
  row.names = FALSE
)

write.csv(
  pp_tab,
  file.path(out_dir, "posterior_predictive_checks.csv"),
  row.names = FALSE
)

cat(
  "bayes_p near 0 or 1 flags a statistic ",
  "the model does not reproduce.\n"
)


## ------------------------------------------------------------------------
## Fitted vs observed by stage
## ------------------------------------------------------------------------

fit_mean <- colMeans(
  
  exp(
    lp_fixed +
      a_draws[
        ,
        match(
          as.character(dat$int_id),
          rat_levels
        )
      ]
  ) *
    
    rep(
      dat$exposure,
      each = S
    )
)


stage_tab <- do.call(
  rbind,
  lapply(
    sort(unique(dat$stage)),
    function(s) {
      
      r <- dat$stage == s
      
      data.frame(
        stage = s,
        
        n = sum(r),
        
        observed_mean =
          mean(dat$y[r]),
        
        fitted_mean =
          mean(fit_mean[r]),
        
        observed_var =
          stats::var(dat$y[r]),
        
        predictive_var =
          mean(
            apply(
              yrep_store[
                ,
                r,
                drop = FALSE
              ],
              1L,
              stats::var
            )
          )
      )
    }
  )
)


stage_tab_print <- stage_tab

num_cols <- vapply(
  stage_tab_print,
  is.numeric,
  logical(1)
)

stage_tab_print[num_cols] <- lapply(
  stage_tab_print[num_cols],
  round,
  digits = 3
)

print(
  stage_tab_print,
  row.names = FALSE
)

write.csv(
  stage_tab,
  file.path(out_dir, "posterior_predictive_checks_by_stage.csv"),
  row.names = FALSE
)


## ------------------------------------------------------------------------
## Posterior predictive plot
## ------------------------------------------------------------------------

pp_plot_df <- rbind(
  
  data.frame(
    count = dat$y,
    distribution = "Observed"
  ),
  
  data.frame(
    count = as.vector(yrep_store),
    distribution = "Posterior predictive"
  )
)


p_pp <- ggplot(
  pp_plot_df,
  aes(x = count)
) +
  geom_histogram(
    aes(
      y = after_stat(density)
    ),
    bins = 40,
    boundary = 0
  ) +
  facet_wrap(
    ~ distribution,
    scales = "free_y"
  ) +
  labs(
    x = "Daily seizure count",
    y = "Density"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

print(p_pp)
ggsave(file.path(out_dir, "posterior_predictive_distribution.pdf"),
       p_pp, width = 7.5, height = 4.8)


## ========================================================================
## 8. PREDICTIVE COUNTS FOR NEW DAYS
## ========================================================================

if (
  file.exists(
    file.path(
      "src",
      "getTheta.cpp"
    )
  ) &&
  requireNamespace(
    "Rcpp",
    quietly = TRUE
  )
) {
  
  if (
    !exists(
      "getTheta_cpp",
      mode = "function"
    )
  ) {
    
    Rcpp::sourceCpp(
      file.path(
        "src",
        "getTheta.cpp"
      )
    )
  }
  
  
  scen_days <- round(
    seq(
      min(day_grid),
      max(day_grid),
      length.out = 5
    )
  )
  
  
  scen <- expand.grid(
    day = scen_days,
    stage = c(1L, 2L)
  )
  
  
  Bs <- rbind(
    
    predict(
      basis,
      t_of_day(
        scen_days,
        1L
      )
    ),
    
    predict(
      basis,
      t_of_day(
        scen_days,
        2L
      )
    )
  )
  
  
  eta_scen <- matrix(
    NA_real_,
    S,
    nrow(scen)
  )
  
  
  for (k in seq_len(nrow(scen))) {
    
    if (scen$stage[k] == 1L) {
      
      eta_scen[, k] <-
        beta_draws[
          ,
          idx_g1,
          drop = FALSE
        ] %*% Bs[k, ] +
        beta_draws[, idx_int]
      
    } else {
      
      eta_scen[, k] <-
        beta_draws[
          ,
          idx_g2,
          drop = FALSE
        ] %*% Bs[k, ] +
        beta_draws[, idx_int] +
        beta_draws[, idx_stage]
    }
  }
  
  
  # One day exposure, typical rat
  mu_scen <- exp(
    eta_scen
  )
  
  
  # Keep means inside support
  mu_scen <- pmin(
    mu_scen,
    max(support) - 0.5
  )
  
  
  pred <- matrix(
    NA_real_,
    length(pp_idx),
    nrow(scen)
  )
  
  
  for (j in seq_along(pp_idx)) {
    
    s <- pp_idx[j]
    
    th <- getTheta_cpp(
      spt = as.numeric(support),
      f0 = w_draws[s, ],
      mu = mu_scen[s, ],
      thetaStart = NULL,
      logit = FALSE
    )$theta
    
    
    pred[j, ] <- draw_counts(
      tilted_pmf(
        th,
        w_draws[s, ]
      ),
      support
    )
  }
  
  
  scen_tab <- data.frame(
    
    scen,
    
    expected =
      colMeans(mu_scen),
    
    pred_median =
      apply(
        pred,
        2L,
        stats::median
      ),
    
    pred_lower =
      apply(
        pred,
        2L,
        stats::quantile,
        0.025
      ),
    
    pred_upper =
      apply(
        pred,
        2L,
        stats::quantile,
        0.975
      ),
    
    `P(zero)` =
      colMeans(pred == 0),
    
    check.names = FALSE
  )
  
  
  cat(
    "\nPredictive daily counts for a typical rat, ",
    "one full day of exposure:\n"
  )
  
  
  # Round numeric columns only
  scen_tab_print <- scen_tab
  
  num_cols <- vapply(
    scen_tab_print,
    is.numeric,
    logical(1)
  )
  
  scen_tab_print[num_cols] <- lapply(
    scen_tab_print[num_cols],
    round,
    digits = 3
  )
  
  print(
    scen_tab_print,
    row.names = FALSE
  )

  write.csv(
    scen_tab,
    file.path(out_dir, "new_day_predictive_counts.csv"),
    row.names = FALSE
  )
  
} else {
  
  cat(
    "\nSkipped new-day prediction ",
    "(needs Rcpp and src/getTheta.cpp).\n"
  )
}


## ========================================================================
## DONE
## ========================================================================

cat(
  "\nDone. Figures and result tables were saved to: ",
  out_dir, "\n"
)
