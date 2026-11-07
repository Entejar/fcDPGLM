#!/usr/bin/env Rscript
# ============================================================================
# Fit fcDPGLM to the seizure count data.
#
# Model: y_i | a_r ~ exponential tilt of a gamma-CRM reference distribution,
# with mean e_i * exp(x_i' beta + a_r), where x_i holds an intercept, a stage
# indicator and stage-specific natural cubic splines in day (3 df each), e_i is
# the exposure in days and a_r ~ N(0, sigma_a^2) is a rat random intercept.
# Priors:
#   intercept and stage effect ~ N(0, 10^2)
#   stage-specific spline coefficients (orthonormalized basis per stage)
#     ~ N(0, tau_d^2), tau_d ~ half-Cauchy(0, 1)
#   sigma_a ~ half-Cauchy(0, 1)
#   alpha ~ Gamma(1, 0.1) (shape, rate)
#   G0 = NB(mean = overall exposure-adjusted rate, size from a preliminary
#     NB GLMM), on the counts 0, ..., K with K = max(y) + 30
# Sampler: Fisher-preconditioned MALA for (beta, a), CRM weights in 4 random
# blocks, an exact Gibbs draw of the total CRM mass and slice sampling of alpha.
# Several chains are run in parallel from dispersed starting values.
#
# Usage (from the repository root):
#   Rscript scripts/01_fit_fcdpglm.R chains=4 iter=30000 adapt=5000 thin=5
# Output: results/fcdpglm_fit.rds
# ============================================================================

suppressPackageStartupMessages({
  library(splines); library(glmmTMB); library(coda)
})

source("R/fcdpglm.R")
# Compile the C++ code once, before the chains are forked (forked workers
# share tempdir(), so compiling inside them races on the same build folder).
rat_dpglm_ensure_getTheta()
rat_dpglm_ensure_kernels()

args <- list(data_path = "data/synthetic_seizure_counts.csv",
             output = "results/fcdpglm_fit.rds",
             chains = 4L, iter = 30000L, adapt = 5000L, thin = 5L, seed = 20261015L,
             mass_move = 1L)
for (a in unlist(strsplit(commandArgs(trailingOnly = TRUE), "[[:space:]]+"))) {
  if (!nzchar(a)) next
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L || !kv[1] %in% names(args)) {
    stop("Unrecognized argument '", a, "'. Valid names: ",
         paste(names(args), collapse = ", "), call. = FALSE)
  }
  args[[kv[1]]] <- if (kv[1] %in% c("chains", "iter", "adapt", "thin", "seed", "mass_move"))
    as.integer(kv[2]) else kv[2]
}
stopifnot(args$adapt < args$iter, args$chains >= 1L, args$thin >= 1L)
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)


# -------------------------------------------------------------------------
# Data: one row per rat-day with columns rat, stage (1 or 2), day (1, 2, ...
# within stage), exposure (observed time in days) and y (seizure count).
# Time within stage is standardized separately for each stage.
# -------------------------------------------------------------------------

df_day <- read.csv(args$data_path)
stopifnot(all(c("rat", "stage", "day", "exposure", "y") %in% names(df_day)),
          all(df_day$exposure > 0), all(df_day$stage %in% 1:2))
df_day$int_id <- as.integer(df_day$rat)
df_day$y <- as.integer(df_day$y)
df_day$t <- ave(df_day$day, df_day$stage, FUN = function(d) as.numeric(scale(d)))
df_day$stage2 <- as.integer(df_day$stage == 2L)
df_day$rat <- factor(df_day$int_id)

basis <- ns(df_day$t, df = 3, intercept = FALSE)
X <- cbind(
  intercept = 1,
  stage2 = df_day$stage2,
  basis * (df_day$stage == 1L),
  basis * (df_day$stage == 2L)
)
colnames(X) <- c("intercept", "stage2",
                 paste0("stage1_spline", seq_len(ncol(basis))),
                 paste0("stage2_spline", seq_len(ncol(basis))))
p <- ncol(X)
idx_g1 <- grep("^stage1_spline", colnames(X)); idx_g2 <- grep("^stage2_spline", colnames(X))
cat("Rat-days:", nrow(df_day), "| rats:", nlevels(df_day$rat),
    "| max count:", max(df_day$y), "\n")


# -------------------------------------------------------------------------
# Preliminary NB GLMM (same parameterization as X): starting values, the
# G0 size, and the initial sigma_a.  Not used as a prior.
# -------------------------------------------------------------------------

design_df <- as.data.frame(X)
names(design_df) <- paste0("X_", make.names(colnames(X)))
model_df <- data.frame(y = df_day$y, exposure = df_day$exposure, rat = df_day$rat,
                       design_df, check.names = FALSE)
eb_formula <- as.formula(paste("y ~ 0 +", paste(names(design_df), collapse = " + "),
                               "+ offset(log(exposure)) + (1 | rat)"))
fit_eb_nb <- glmmTMB(eb_formula, data = model_df, family = nbinom2(link = "log"))
beta_eb <- setNames(as.numeric(fixef(fit_eb_nb)$cond), colnames(X))
V_eb    <- as.matrix(vcov(fit_eb_nb)$cond)
a_eb    <- ranef(fit_eb_nb)$cond$rat[as.character(unique(df_day$int_id)), 1]
sigma_a_eb <- sqrt(VarCorr(fit_eb_nb)$cond$rat[1, 1])
g0_size <- sigma(fit_eb_nb)
g0_mu   <- sum(df_day$y) / sum(df_day$exposure)
cat("NB GLMM: sigma_a =", round(sigma_a_eb, 4), "| NB size =", round(g0_size, 4),
    "| G0 mean =", round(g0_mu, 3), "\n")


# -------------------------------------------------------------------------
# Orthonormalized spline basis per stage:
# X_fit = X T, beta = T beta_fit.  Draws are mapped back to the original basis.
# -------------------------------------------------------------------------

T_hier <- diag(p)
for (s in 1:2) {
  g <- if (s == 1) idx_g1 else idx_g2
  Bd <- X[df_day$stage == s, g, drop = FALSE]
  T_hier[g, g] <- solve(qr.R(qr(Bd))) * sqrt(nrow(Bd))
}
X_fit <- X %*% T_hier


# -------------------------------------------------------------------------
# Chains: chain 1 starts at the NB GLMM estimates; the others at dispersed
# values (coefficients drawn from N(beta_eb, 4 V_eb), rat effects jittered).
# -------------------------------------------------------------------------

run_chain <- function(ch) {
  set.seed(args$seed + ch)
  b0 <- if (ch == 1L) beta_eb else MASS::mvrnorm(1, beta_eb, 4 * (V_eb + t(V_eb)) / 2)
  a0 <- if (ch == 1L) a_eb else a_eb + rnorm(length(a_eb), 0, 0.5 * sigma_a_eb)
  t0 <- proc.time()[3]
  fit <- fit_rat_dpglm(
    y = df_day$y, X = X_fit, exposure = df_day$exposure, rat_id = df_day$int_id,
    iter = args$iter, adapt = args$adapt,
    g0 = "nbinom", g0_mu = g0_mu, g0_size = g0_size,
    max_count = max(df_day$y) + 30L, warmup = 25L,
    alpha = 1, alpha_prior = c(1, 0.1),
    sigma_a = sigma_a_eb, sigma_a_prior = "halfcauchy", sigma_a_scale = 1,
    mu_beta = rep(0, p), sigma_beta = rep(10, p),
    beta_groups = list(idx_g1, idx_g2), beta_group_scale = 1,
    initial_beta = drop(solve(T_hier, b0)), initial_a = a0,
    MCMC_update = "fisher_mala", crm_blocks = 4L,
    theta_warm_start = TRUE, use_cpp_kernels = TRUE,
    crm_mass_move = args$mass_move == 1L,   # exact Gibbs draw of the CRM total mass
    alpha_update = if (args$mass_move == 1L) "conditional" else "joint",
    report_progress = (ch == 1L))
  fit$beta <- fit$beta %*% t(T_hier)          # original basis
  colnames(fit$beta) <- colnames(X)
  # keep the post-adaptation draws only, thinned, to limit memory
  keep <- seq.int(args$adapt + 1L, args$iter, by = args$thin)
  list(chain = ch, seconds = proc.time()[3] - t0,
       beta = fit$beta[keep, , drop = FALSE], a = fit$a[keep, , drop = FALSE],
       theta = fit$theta[keep, , drop = FALSE], u = fit$u[keep, , drop = FALSE],
       weights = fit$weights[keep, , drop = FALSE], loglik = fit$loglik[keep],
       crm_log_r = fit$crm_log_r[keep],
       sigma_a = fit$sigma_a_samples[keep], alpha = fit$alpha_samples[keep],
       tau = fit$tau_samples[keep, , drop = FALSE],
       acceptance = fit$acceptance, tuning = fit$tuning, support = fit$support,
       shape = fit$shape, g0 = fit$g0)
}

cat("\nRunning", args$chains, "chains of", args$iter, "iterations (adaptation",
    args$adapt, ", thin", args$thin, ")\n")
chains <- parallel::mclapply(seq_len(args$chains), function(ch)
  tryCatch(run_chain(ch), error = function(e) e),
  mc.cores = min(args$chains, max(1L, parallel::detectCores() - 1L)),
  mc.preschedule = FALSE)
failed <- vapply(chains, inherits, logical(1), "error")
if (any(failed)) {
  for (i in which(failed)) cat("Chain", i, "failed:", conditionMessage(chains[[i]]), "\n")
  chains <- chains[!failed]
  if (!length(chains)) stop("All chains failed.", call. = FALSE)
}


# -------------------------------------------------------------------------
# Convergence diagnostics across chains (R-hat, ESS) for key quantities
# -------------------------------------------------------------------------

functionals <- function(ch) {
  b <- ch$beta
  cbind(b, sigma_a = ch$sigma_a, log_alpha = log(ch$alpha), loglik = ch$loglik,
        crm_mass = log(rowSums(ch$weights)))
}
mcl <- coda::mcmc.list(lapply(chains, function(ch) coda::mcmc(functionals(ch))))
rhat <- if (length(chains) > 1L) coda::gelman.diag(mcl, autoburnin = FALSE,
                                                   multivariate = FALSE)$psrf[, 1] else NA
ess  <- coda::effectiveSize(mcl)
diag_tab <- data.frame(parameter = colnames(functionals(chains[[1]])),
                       rhat = round(rhat, 3), ess = round(ess))
cat("\nConvergence diagnostics (R-hat < 1.01 and ESS > 400 are good):\n")
print(diag_tab, row.names = FALSE)
cat("\nPost-adaptation acceptance rates by chain:\n")
print(t(sapply(chains, function(ch) unlist(ch$acceptance$post_adaptation[c("joint", "crm", "alpha")]))))
cat("Run time (minutes) by chain:", round(sapply(chains, `[[`, "seconds") / 60, 1), "\n")


# -------------------------------------------------------------------------
# Pooled fit in the layout read by scripts/03_posterior_summaries.R:
# all retained draws, adaptation already removed (tuning$adapt_until = 0).
# -------------------------------------------------------------------------

pool <- function(name) do.call(rbind, lapply(chains, `[[`, name))
pooled <- list(
  model = "fcDPGLM (pooled chains)",
  beta = pool("beta"), a = pool("a"), theta = pool("theta"), u = pool("u"),
  weights = pool("weights"), loglik = unlist(lapply(chains, `[[`, "loglik")),
  crm_log_r = unlist(lapply(chains, `[[`, "crm_log_r")),
  sigma_a_samples = unlist(lapply(chains, `[[`, "sigma_a")),
  alpha_samples = unlist(lapply(chains, `[[`, "alpha")),
  tau_samples = pool("tau"),
  sigma_a = median(unlist(lapply(chains, `[[`, "sigma_a"))),
  support = chains[[1]]$support, shape = chains[[1]]$shape, g0 = chains[[1]]$g0,
  chain_id = rep(vapply(chains, `[[`, integer(1), "chain"),
                 vapply(chains, function(ch) nrow(ch$beta), integer(1))),
  # modifyList (not c()) so adapt_until is replaced: the pooled draws are all
  # post-adaptation and none should be discarded by the analysis script.
  tuning = modifyList(chains[[1]]$tuning, list(adapt_until = 0L, thin = args$thin,
                                               chains = length(chains))),
  acceptance = {
    acc_mean <- function(f) mean(vapply(chains, function(ch) {
      v <- f(ch$acceptance); if (is.null(v) || !length(v)) NA_real_ else as.numeric(v) }, numeric(1)))
    list(beta = acc_mean(function(x) x$beta), random_effects = acc_mean(function(x) x$random_effects),
         joint = acc_mean(function(x) x$joint), crm = acc_mean(function(x) x$crm),
         post_adaptation = list(
           joint = acc_mean(function(x) x$post_adaptation$joint),
           crm = acc_mean(function(x) x$post_adaptation$crm),
           alpha = acc_mean(function(x) x$post_adaptation$alpha)))
  }
)

saveRDS(list(
  data = df_day, design_matrix = X, glmmTMB = fit_eb_nb, rat_dpglm = pooled,
  chains_summary = list(rhat_ess = diag_tab,
                        seconds = sapply(chains, `[[`, "seconds")),
  T_hier = T_hier,
  settings = c(args, list(
    prior = paste("alpha ~ Gamma(1, 0.1); sigma_a ~ half-Cauchy(0, 1);",
                  "intercept, stage ~ N(0, 10^2); orthonormalized stage-specific",
                  "spline coefficients ~ N(0, tau_d^2), tau_d ~ half-Cauchy(0, 1)"),
    g0 = sprintf("NB(mean %.4f, size %.4f), K = %d", g0_mu, g0_size, max(df_day$y) + 30L),
    sampler = "Fisher-preconditioned MALA for (beta, a); CRM in 4 blocks; exact total-mass draw and slice update of alpha"))
), args$output)

cat("\nSaved", args$output, "with", nrow(pooled$beta), "pooled draws from",
    length(chains), "chains.\n")
cat("alpha: posterior median", round(median(pooled$alpha_samples), 1),
    "| sigma_a: posterior median", round(median(pooled$sigma_a_samples), 3), "\n")
