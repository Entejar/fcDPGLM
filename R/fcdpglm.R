# fcDPGLM: a semiparametric Bayesian generalized linear mixed model with a
# nonparametric count distribution, for repeated count data.
#
# The reference distribution is a normalized gamma completely random measure on
# the counts 0, ..., K, with K well above the largest observed count, i.e.
# independent gamma weights w_k ~ Gamma(alpha * G0(k), 1).  After the
# auxiliary-u augmentation the full conditional of the weights is conjugate,
#
#   w_k | u, theta, y ~ Gamma(alpha * G0(k) + n*_k, rate = 1 + psi(k)),
#
# and the normalising constant of this proposal is available in closed form,
# so the Metropolis-Hastings ratio for the CRM update is exact.  An MH step is
# needed because theta is a deterministic function of the weights through the
# mean constraint b'(theta_i) = mu_i.

# Compile src/getTheta.cpp unless a working getTheta_cpp() is already loaded
# (an Rcpp function restored from a saved workspace has an invalid pointer, so
# it is called once to check).
rat_dpglm_ensure_getTheta <- function() {
  works <- exists("getTheta_cpp", mode = "function") && isTRUE(tryCatch({
    probe <- getTheta_cpp(spt = c(0, 1, 2), f0 = c(1, 1, 1), mu = 1,
                          thetaStart = NULL, logit = FALSE)$theta
    length(probe) == 1L && is.finite(probe)
  }, error = function(error) FALSE))
  if (!works) {
    cpp_file <- file.path("src", "getTheta.cpp")
    if (!file.exists(cpp_file)) {
      stop("Cannot find ", cpp_file,
           ". Set the working directory to the repository root.", call. = FALSE)
    }
    Rcpp::sourceCpp(cpp_file)
  }
  invisible(NULL)
}

# Options set by fit_rat_dpglm() for the duration of a fit (restored on exit):
#   cpp: use the C++ kernels in src/fcdpglm_kernels.cpp for b(theta), log psi
#        and the tilted variances (same quantities, faster).
.rat_dpglm_opts <- new.env()
.rat_dpglm_opts$cpp <- FALSE

rat_dpglm_ensure_kernels <- function() {
  works <- exists("rat_btheta_cpp", mode = "function") && isTRUE(tryCatch(
    is.finite(rat_btheta_cpp(0, c(0, 1), c(0, 0))), error = function(error) FALSE))
  if (!works) {
    cpp_file <- file.path("src", "fcdpglm_kernels.cpp")
    if (!file.exists(cpp_file)) {
      stop("Cannot find ", cpp_file,
           ". Set the working directory to the repository root.", call. = FALSE)
    }
    Rcpp::sourceCpp(cpp_file)
  }
  invisible(NULL)
}

# Describe why the log posterior is not finite at a given state.
rat_dpglm_diagnose_state <- function(beta, random_effects, weights, support, y,
                                     X, exposure, rat_index, mu_beta,
                                     sigma_beta, sigma_a) {
  lines <- character(0)
  mean_values <- exposure * exp(drop(X %*% beta) + random_effects[rat_index])
  lines <- c(lines, sprintf(
    "fitted means run %.4g to %.4g; the support is 0 to %d",
    min(mean_values), max(mean_values), max(support)))
  lines <- c(lines, sprintf(
    "weights: %d non-finite, %d not positive, range %.4g to %.4g",
    sum(!is.finite(weights)), sum(weights <= 0), min(weights), max(weights)))
  state <- tryCatch(
    rat_dpglm_state(beta, random_effects, weights, support, X, exposure, rat_index),
    error = function(error) conditionMessage(error)
  )
  if (is.character(state)) {
    lines <- c(lines, paste("the theta solve failed with:", state))
  } else {
    lines <- c(lines, sprintf(
      "theta: %d non-finite, range %.4g to %.4g",
      sum(!is.finite(state$theta)), min(state$theta), max(state$theta)))
    lines <- c(lines, sprintf(
      "b(theta): %d non-finite, range %.4g to %.4g",
      sum(!is.finite(state$btheta)), min(state$btheta), max(state$btheta)))
    lines <- c(lines, sprintf(
      "smallest weight at an observed count: %.4g", min(weights[y + 1L])))
    lines <- c(lines, paste(
      "log likelihood:", format(rat_dpglm_loglik(y, state, weights))))
  }
  lines <- c(lines, paste(
    "beta prior:", format(sum(dnorm(beta, mu_beta, sigma_beta, log = TRUE))),
    "| random-effect prior:",
    format(sum(dnorm(random_effects, 0, sigma_a, log = TRUE)))))
  paste(lines, collapse = "\n  ")
}

# log(1 + exp(log_x)) without overflowing for large log_x.
rat_dpglm_log1p_exp <- function(log_x) {
  out <- log_x
  small <- !is.na(log_x) & log_x < 30
  out[small] <- log1p(exp(log_x[small]))
  out
}

# b_i(theta_i; w) = log sum_k w_k exp(theta_i * k), one value per observation.
rat_dpglm_btheta <- function(theta, atoms, weights) {
  if (.rat_dpglm_opts$cpp) return(rat_btheta_cpp(theta, atoms, log(weights)))
  terms <- sweep(outer(theta, atoms), 2L, log(weights), "+")
  max_terms <- apply(terms, 1L, max)
  max_terms + log(rowSums(exp(terms - max_terms)))
}

# log psi(k) = log sum_i u_i exp(theta_i * k), one value per support point.
rat_dpglm_log_psi <- function(log_u, theta, support) {
  if (.rat_dpglm_opts$cpp) return(rat_log_psi_cpp(log_u, theta, support))
  terms <- outer(theta, support) + log_u
  max_terms <- apply(terms, 2L, max)
  max_terms + log(colSums(exp(sweep(terms, 2L, max_terms, "-"))))
}

# theta_start: optional warm start for the Newton solve of b'(theta) = mean
# (e.g. the current theta); NULL starts every solve at 0.
rat_dpglm_state <- function(beta, random_effects, weights, support, X, exposure,
                            rat_index, theta_start = NULL) {
  mean_values <- exposure * exp(drop(X %*% beta) + random_effects[rat_index])
  if (!is.null(theta_start) &&
      (length(theta_start) != length(mean_values) || any(!is.finite(theta_start)))) {
    theta_start <- NULL
  }
  theta <- getTheta_cpp(
    spt = support,
    f0 = weights,
    mu = mean_values,
    thetaStart = theta_start,
    logit = FALSE
  )$theta
  list(
    mean = mean_values,
    theta = theta,
    btheta = rat_dpglm_btheta(theta, support, weights)
  )
}

# Log-likelihood contribution of one rat's rows.  The log w_y terms do not
# involve the rat effect, so they are dropped: this is used only inside the
# random-effect Metropolis ratio, where they cancel.
rat_dpglm_block_loglik <- function(a_value, beta, rows, weights, support, X,
                                   exposure, y) {
  mean_values <- exposure[rows] *
    exp(drop(X[rows, , drop = FALSE] %*% beta) + a_value)
  block <- tryCatch({
    theta <- getTheta_cpp(
      spt = support, f0 = weights, mu = mean_values,
      thetaStart = NULL, logit = FALSE
    )$theta
    list(theta = theta, btheta = rat_dpglm_btheta(theta, support, weights))
  }, error = function(error) NULL)
  if (is.null(block) || any(!is.finite(block$theta)) ||
      any(!is.finite(block$btheta))) {
    return(-Inf)
  }
  sum(block$theta * y[rows] - block$btheta)
}

rat_dpglm_loglik <- function(y, state, weights) {
  weight_at_y <- weights[y + 1L]
  if (any(!is.finite(weight_at_y)) || any(weight_at_y <= 0)) {
    return(-Inf)
  }
  sum(state$theta * y - state$btheta + log(weight_at_y))
}

rat_dpglm_logposterior <- function(beta, random_effects, weights, support, y, X,
                                   exposure, rat_index, mu_beta, sigma_beta,
                                   sigma_a, theta_start = NULL) {
  if (any(!is.finite(beta)) || any(!is.finite(random_effects)) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    return(-Inf)
  }
  state <- tryCatch(
    rat_dpglm_state(
      beta, random_effects, weights, support, X, exposure, rat_index, theta_start
    ),
    error = function(error) NULL
  )
  if (is.null(state) || any(!is.finite(state$theta))) {
    return(-Inf)
  }
  loglik <- rat_dpglm_loglik(y, state, weights)
  log_beta_prior <- sum(dnorm(beta, mu_beta, sigma_beta, log = TRUE))
  log_a_prior <- sum(dnorm(random_effects, 0, sigma_a, log = TRUE))
  loglik + log_beta_prior + log_a_prior
}

# Draw from the conjugate full conditional of the weights given (u, theta).
rat_dpglm_sample_weights <- function(shape, log_psi, floor_value) {
  scale <- exp(-rat_dpglm_log1p_exp(log_psi))
  pmax(stats::rgamma(length(shape), shape = shape, scale = scale), floor_value)
}

# Exact log Metropolis-Hastings ratio for the CRM update.
#
#   log R = sum_i (theta*_i - theta_i) y_i
#           + sum_i u_i [T_i(w*) + T_i(w) - T*_i(w*) - T*_i(w)]
#           + sum_k shape_k [log(1 + psi*(k)) - log(1 + psi(k))]
#
# with T_i(w) = exp(b_i(theta_i; w)).  The first line is the likelihood, the
# second comes from the exponential part of the gamma weights (it is the
# sum_k (w*_k + w_k) (psi(k) - psi*(k)) term written observation-wise), and the
# third is the ratio of the proposal normalising constants.  When theta* equals
# theta the ratio is exactly zero, i.e. the step reduces to a Gibbs draw.
rat_dpglm_crm_log_ratio <- function(y, log_u, shape, theta, theta_star,
                                    support, weights, proposed_weights,
                                    log_psi, log_psi_star, btheta_current,
                                    btheta_star_proposed) {
  b_current_proposed <- rat_dpglm_btheta(theta, support, proposed_weights)
  b_star_current <- rat_dpglm_btheta(theta_star, support, weights)
  
  data_term <- sum((theta_star - theta) * y)
  mass_term <- sum(
    exp(log_u + b_current_proposed) + exp(log_u + btheta_current) -
      exp(log_u + btheta_star_proposed) - exp(log_u + b_star_current)
  )
  laplace_term <- sum(
    shape * (rat_dpglm_log1p_exp(log_psi_star) - rat_dpglm_log1p_exp(log_psi))
  )
  
  total_term <- data_term + mass_term + laplace_term
  attr(total_term, "data_term") <- data_term
  attr(total_term, "mass_term") <- mass_term
  attr(total_term, "laplace_term") <- laplace_term
  total_term
}

# Exact log of draws from Gamma(shape, rate = 1 + psi(k)), accurate when the
# shape is tiny (where rgamma underflows to 0): G = G1 * U^(1 / shape) with
# G1 ~ Gamma(shape + 1, 1), U ~ Uniform(0, 1).
rat_dpglm_sample_log_weights <- function(shape, log_psi) {
  log(stats::rgamma(length(shape), shape = shape + 1, rate = 1)) +
    log(stats::runif(length(shape))) / shape - rat_dpglm_log1p_exp(log_psi)
}

# Exact log MH ratio for a CRM update that changes only the weights in `block`
# (a subset of support indices), with u held fixed.  With a_k = w*_k + w_k for
# k in the block and a_k = w_k elsewhere, the exponential terms collapse to
#   sum_i u_i [T_i(theta; a) - T_i(theta*; a)],
# and only block points enter the Laplace term.  With the whole support as the
# block this equals rat_dpglm_crm_log_ratio().
rat_dpglm_crm_block_log_ratio <- function(y, log_u, shape, theta, theta_star,
                                          support, weights, proposed_weights,
                                          log_psi, log_psi_star, block) {
  mixed <- weights
  mixed[block] <- proposed_weights[block] + weights[block]
  data_term <- sum((theta_star - theta) * y)
  mass_term <- sum(exp(log_u + rat_dpglm_btheta(theta, support, mixed)) -
                     exp(log_u + rat_dpglm_btheta(theta_star, support, mixed)))
  laplace_term <- sum(shape[block] * (rat_dpglm_log1p_exp(log_psi_star[block]) -
                                        rat_dpglm_log1p_exp(log_psi[block])))
  total_term <- data_term + mass_term + laplace_term
  attr(total_term, "data_term") <- data_term
  attr(total_term, "mass_term") <- mass_term
  attr(total_term, "laplace_term") <- laplace_term
  total_term
}

# Gradient and Fisher information of the log posterior in xi = (beta, a) for
# the Fisher-preconditioned MALA update, holding the weights fixed.
#   d loglik / d eta_i = (y_i - m_i) m_i / V_i, with V_i = b''(theta_i), because
#   d theta / d m = 1 / V and d m / d eta = m; Fisher weights are m_i^2 / V_i.
rat_dpglm_mala_terms <- function(beta, random_effects, state, weights, support,
                                 y, Z, mu_beta, sigma_beta, sigma_a) {
  m <- state$mean
  v <- if (.rat_dpglm_opts$cpp) {
    rat_tilted_var_cpp(state$theta, support, log(weights), state$btheta)
  } else {
    pr <- exp(outer(state$theta, support) +
                rep(log(weights), each = length(state$theta)) - state$btheta)
    drop(pr %*% support^2) - drop(pr %*% support)^2
  }
  v <- pmax(v, 1e-10)
  score <- (y - m) * m / v
  prior_prec <- c(1 / sigma_beta^2, rep(1 / sigma_a^2, length(random_effects)))
  xi <- c(beta, random_effects)
  grad <- drop(crossprod(Z, score)) - prior_prec * (xi - c(mu_beta, rep(0, length(random_effects))))
  info <- crossprod(Z, Z * (m^2 / v)) + diag(prior_prec, length(xi))
  R <- tryCatch(chol((info + t(info)) / 2), error = function(error) NULL)
  list(grad = grad, R = R)
}

# log density (up to a constant common to both directions) of the MALA proposal
# xi_to ~ N(xi_from + (eps^2 / 2) H^-1 g, eps^2 H^-1), with H = R'R.
rat_dpglm_mala_logq <- function(xi_to, xi_from, terms, eps) {
  drift <- backsolve(terms$R, forwardsolve(t(terms$R), terms$grad))
  r <- terms$R %*% (xi_to - xi_from - 0.5 * eps^2 * drift)
  -0.5 * sum(r^2) / eps^2 + sum(log(diag(terms$R)))
}

fit_rat_dpglm <- function(y, X, exposure, rat_id, iter = 2000L,
                          alpha = 1, g0 = c("nbinom", "poisson"),
                          g0_mu = NULL, g0_size = 10, max_count = NULL,
                          mu_beta = NULL, sigma_beta = NULL, sigma_a = NULL,
                          initial_beta = NULL, initial_a = NULL,
                          initial_weights = NULL,
                          beta_step = 0.05, beta_cov = NULL, a_step = 0.15,
                          joint_cov = NULL, joint_step = NULL,
                          warmup = 25L, adapt = NULL, report_progress = TRUE,
                          MCMC_update = c("rw", "fisher_mala"), mala_step = 1,
                          crm_blocks = 1L,
                          sigma_a_prior = c("fixed", "halfcauchy", "invgamma"),
                          sigma_a_scale = 1, sigma_a_ig = c(1, 1),
                          alpha_prior = NULL, alpha_step = 0.5,
                          beta_groups = NULL, beta_group_scale = 1,
                          theta_warm_start = FALSE, use_cpp_kernels = FALSE,
                          crm_mass_move = FALSE,
                          alpha_update = c("joint", "conditional")) {
  # MCMC_update selects the (beta, a) update:
  #   "rw"          random-walk Metropolis (joint if joint_cov is
  #                 given, otherwise beta block + scalar a_r updates).
  #   "fisher_mala" joint (beta, a) Metropolis-adjusted Langevin proposal
  #                 preconditioned by the Fisher information at the current
  #                 state; mala_step is adapted toward acceptance 0.574.
  # crm_blocks: 1 = update all CRM weights at once; > 1 splits
  # the support at random into that many blocks each iteration and applies the
  # same conjugate-proposal MH step to each block in turn.
  #
  # Optional hyperparameter updates (by default these are held fixed):
  #   sigma_a_prior  "fixed" (sigma_a held at its supplied value), "halfcauchy"
  #                  (sigma_a ~ half-Cauchy(0, sigma_a_scale), Gibbs via the
  #                  inverse-gamma mixture of Makalic and Schmidt 2016), or
  #                  "invgamma" (sigma_a^2 ~ IG(sigma_a_ig[1], sigma_a_ig[2]),
  #                  rate parameterisation).  sigma_a is then the initial value.
  #   alpha_prior    NULL (alpha fixed) or c(shape, rate) of a gamma prior; alpha
  #                  is updated jointly with all CRM weights by an exact MH step
  #                  (random walk on log alpha with sd alpha_step, adapted toward
  #                  acceptance 0.25; weights from the conjugate proposal).
  #   beta_groups    NULL, or a list of integer column groups of X; the
  #                  coefficients in group g get N(0, tau_g^2) priors with
  #                  tau_g ~ half-Cauchy(0, beta_group_scale), Gibbs-updated.
  #                  Columns outside the groups keep N(mu_beta, sigma_beta^2).
  # alpha_update   how alpha is updated when alpha_prior is given: "joint"
  #                (a joint MH move of alpha and all weights; it mixes
  #                slowly because it is only accepted for tiny changes in alpha)
  #                or "conditional" (slice sampling of log alpha from its exact
  #                conditional given the weights,
  #                  p(alpha | w) propto p(alpha) prod_k Gamma(w_k; alpha G0(k), 1),
  #                using exact log-weights: weights at unobserved counts are far
  #                below the double-precision floor, so they are tracked on the
  #                log scale).  Use with crm_mass_move = TRUE.
  # crm_mass_move  TRUE adds an exact Gibbs draw of the total CRM mass
  #                M = sum(w) each iteration: w -> c w and u -> u / c with
  #                c = M*/M, M* ~ Gamma(alpha, 1).  theta and the likelihood are
  #                unchanged by this rescaling, and given the normalized weights
  #                and u * M the conditional of M is exactly Gamma(alpha, 1).
  #                This improves the mixing of M and of alpha.
  # Speed options (same target, same quantities):
  #   theta_warm_start  start each Newton solve for theta at the current theta.
  #   use_cpp_kernels   use src/fcdpglm_kernels.cpp for b(theta), log psi and
  #                     the tilted variances.
  MCMC_update <- match.arg(MCMC_update)
  sigma_a_prior <- match.arg(sigma_a_prior)
  alpha_update <- match.arg(alpha_update)
  if (!is.null(alpha_prior) && (length(alpha_prior) != 2L || any(!is.finite(alpha_prior)) ||
                                any(alpha_prior <= 0))) {
    stop("alpha_prior must be NULL or c(shape, rate) with positive entries.", call. = FALSE)
  }
  old_cpp <- .rat_dpglm_opts$cpp
  on.exit(.rat_dpglm_opts$cpp <- old_cpp, add = TRUE)
  if (isTRUE(use_cpp_kernels)) {
    rat_dpglm_ensure_kernels()
    .rat_dpglm_opts$cpp <- TRUE
  } else {
    .rat_dpglm_opts$cpp <- FALSE
  }
  theta_warm_start <- isTRUE(theta_warm_start)
  crm_blocks <- as.integer(crm_blocks)
  if (length(crm_blocks) != 1L || is.na(crm_blocks) || crm_blocks < 1L) {
    stop("crm_blocks must be a positive integer.", call. = FALSE)
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Package 'Rcpp' is required.", call. = FALSE)
  }
  rat_dpglm_ensure_getTheta()
  
  y <- as.integer(y)
  X <- as.matrix(X)
  exposure <- as.numeric(exposure)
  rat_id <- as.character(rat_id)
  if (length(y) != nrow(X) || length(y) != length(exposure) ||
      length(y) != length(rat_id)) {
    stop("y, exposure, rat_id, and nrow(X) must have matching lengths.", call. = FALSE)
  }
  if (anyNA(y) || any(y < 0) || anyNA(X) || anyNA(exposure) ||
      any(!is.finite(exposure)) || any(exposure <= 0)) {
    stop("y must be nonnegative counts; X and exposure must be finite; exposure must be positive.", call. = FALSE)
  }
  if (length(iter) != 1L || !is.finite(iter) || iter < 2) {
    stop("iter must be at least 2.", call. = FALSE)
  }
  iter <- as.integer(iter)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be positive.", call. = FALSE)
  }
  g0 <- match.arg(g0)
  if (length(g0_size) != 1L || !is.finite(g0_size) || g0_size <= 0) {
    stop("g0_size must be positive.", call. = FALSE)
  }
  if (length(warmup) != 1L || !is.finite(warmup) || warmup < 0) {
    stop("warmup must be a nonnegative number of sweeps.", call. = FALSE)
  }
  warmup <- as.integer(warmup)
  # Step sizes are tuned up to this iteration and held fixed afterwards.
  adapt_until <- if (is.null(adapt)) floor(iter / 2) else as.integer(adapt)
  if (length(adapt_until) != 1L || !is.finite(adapt_until) || adapt_until < 0) {
    stop("adapt must be a nonnegative iteration count.", call. = FALSE)
  }
  
  n <- length(y)
  p <- ncol(X)
  rat_levels <- unique(rat_id)
  rat_index <- match(rat_id, rat_levels)
  n_rats <- length(rat_levels)
  
  if (is.null(g0_mu)) {
    g0_mu <- sum(y) / sum(exposure)
  }
  if (!is.finite(g0_mu) || g0_mu <= 0) {
    stop("g0_mu must be positive.", call. = FALSE)
  }
  if (is.null(max_count)) {
    max_count <- max(max(y), ceiling(4 * g0_mu + 10))
  }
  if (length(max_count) != 1L || !is.finite(max_count) || max_count < max(y)) {
    stop("max_count must be at least max(y).", call. = FALSE)
  }
  max_count <- as.integer(max_count)
  support <- 0:max_count
  
  if (g0 == "nbinom") {
    g0_prob <- stats::dnbinom(support, size = g0_size, mu = g0_mu)
  } else {
    g0_prob <- stats::dpois(support, lambda = g0_mu)
  }
  g0_prob <- g0_prob / sum(g0_prob)
  
  if (is.null(mu_beta)) mu_beta <- rep(0, p)
  if (is.null(sigma_beta)) sigma_beta <- rep(2, p)
  if (length(mu_beta) != p || length(sigma_beta) != p ||
      any(!is.finite(mu_beta)) || any(!is.finite(sigma_beta)) ||
      any(sigma_beta <= 0)) {
    stop("mu_beta and sigma_beta must be finite vectors of length ncol(X); sigma_beta must be positive.", call. = FALSE)
  }
  if (is.null(sigma_a)) sigma_a <- 1
  if (length(sigma_a) != 1L || !is.finite(sigma_a) || sigma_a <= 0) {
    stop("sigma_a must be positive.", call. = FALSE)
  }
  if (!is.null(beta_groups)) {
    beta_groups <- lapply(beta_groups, as.integer)
    if (any(vapply(beta_groups, function(g) any(g < 1L | g > p), logical(1))) ||
        anyDuplicated(unlist(beta_groups))) {
      stop("beta_groups must be disjoint sets of column indices of X.", call. = FALSE)
    }
    for (g in beta_groups) {
      mu_beta[g] <- 0
      sigma_beta[g] <- beta_group_scale
    }
  }
  # Auxiliary variables of the half-Cauchy mixture representations.
  xi_sigma_a <- 1
  xi_beta_groups <- rep(1, length(beta_groups))

  # Random-walk proposal for beta: along beta_cov if supplied (for example the
  # GLMM covariance), with beta_step a scale multiplier defaulting to
  # 2.4/sqrt(p); otherwise independent steps of size beta_step.
  beta_chol <- NULL
  if (!is.null(beta_cov)) {
    beta_cov <- as.matrix(beta_cov)
    if (!identical(dim(beta_cov), c(p, p)) || any(!is.finite(beta_cov))) {
      stop("beta_cov must be a finite ncol(X) by ncol(X) matrix.", call. = FALSE)
    }
    beta_chol <- tryCatch(
      chol((beta_cov + t(beta_cov)) / 2),
      error = function(error) NULL
    )
    if (is.null(beta_chol)) {
      warning("beta_cov is not positive definite; using the independent random walk.",
              call. = FALSE)
    } else if (missing(beta_step)) {
      beta_step <- 2.4 / sqrt(p)
    }
  }
  

  # Optional joint random-walk proposal for xi = (beta, a).
  # If joint_cov is supplied, this replaces the separate beta proposal and
  # sequential scalar random-effect proposals. joint_cov should normally be a
  # regularised covariance matrix estimated from a pilot chain.
  joint_chol <- NULL
  joint_dim <- p + n_rats
  if (!is.null(joint_cov)) {
    joint_cov <- as.matrix(joint_cov)
    if (!identical(dim(joint_cov), c(joint_dim, joint_dim)) ||
        any(!is.finite(joint_cov))) {
      stop("joint_cov must be a finite (p + n_rats) by (p + n_rats) matrix.",
           call. = FALSE)
    }
    joint_cov <- (joint_cov + t(joint_cov)) / 2
    joint_chol <- tryCatch(chol(joint_cov), error = function(error) NULL)
    if (is.null(joint_chol)) {
      stop("joint_cov is not positive definite; regularise it before fitting.",
           call. = FALSE)
    }
    if (is.null(joint_step)) joint_step <- 2.4 / sqrt(joint_dim)
    if (length(joint_step) != 1L || !is.finite(joint_step) || joint_step <= 0) {
      stop("joint_step must be a positive finite scalar.", call. = FALSE)
    }
  }

  if (is.null(initial_beta)) {
    beta <- tryCatch(
      stats::coef(stats::glm.fit(X, y, family = stats::poisson(), offset = log(exposure))),
      error = function(error) rep(0, p)
    )
    beta[!is.finite(beta)] <- 0
  } else {
    beta <- as.numeric(initial_beta)
    if (length(beta) != p) stop("initial_beta has the wrong length.", call. = FALSE)
  }
  if (is.null(initial_a)) {
    random_effects <- rep(0, n_rats)
  } else {
    random_effects <- as.numeric(initial_a)
    if (length(random_effects) != n_rats) stop("initial_a has the wrong length.", call. = FALSE)
  }
  
  # Check that the starting values give fitted means on the scale of the data
  # (coefficients from a different parameterisation of X would not).
  initial_mean <- exposure * exp(drop(X %*% beta) + random_effects[rat_index])
  mean_ratio <- mean(initial_mean) / mean(y)
  if (!is.finite(mean_ratio) || mean_ratio < 1 / 3 || mean_ratio > 3) {
    warning(sprintf(
      paste("initial_beta implies a mean count of %.3g against an observed mean of %.3g.",
            "Check that initial_beta and mu_beta are stated in the same",
            "parameterisation as the columns of X."),
      mean(initial_mean), mean(y)
    ), call. = FALSE)
  }

  zstar <- sort(unique(y))
  nstar <- tabulate(match(y, zstar), nbins = length(zstar))
  # Multiplicity of every support point, and the conjugate gamma shapes.
  nstar_support <- tabulate(y + 1L, nbins = length(support))
  shape <- alpha * g0_prob + nstar_support
  # Every support point keeps a strictly positive weight so that the mean
  # inversion always has the full support available; the floor is at the
  # smallest normal double and is numerically irrelevant.
  weight_floor <- .Machine$double.xmin
  
  # Initial weights: a draw from the conjugate conditional with psi = 0 has
  # total mass about n + alpha, far above the posterior, so the u and w
  # conditionals are cycled `warmup` times (without accept/reject) at the
  # initial (beta, a), keeping the state with the highest log posterior.
  # This affects only the starting values.
  weights <- pmax(stats::rgamma(length(support), shape = shape, rate = 1), weight_floor)
  # Supplied initial weights are used as given, without warm-up sweeps.
  if (!is.null(initial_weights)) {
    initial_weights <- as.numeric(initial_weights)
    if (length(initial_weights) != length(support) ||
        any(!is.finite(initial_weights)) || any(initial_weights < 0)) {
      stop("initial_weights must be max_count + 1 nonnegative finite values.", call. = FALSE)
    }
    weights <- pmax(initial_weights, weight_floor)
    warmup <- 0L
  }
  best_weights <- weights
  best_lp <- rat_dpglm_logposterior(
    beta, random_effects, weights, support, y, X, exposure, rat_index,
    mu_beta, sigma_beta, sigma_a
  )
  if (!is.finite(best_lp)) {
    stop("Initial state has non-finite posterior.\n  ",
         rat_dpglm_diagnose_state(
           beta, random_effects, weights, support, y, X, exposure, rat_index,
           mu_beta, sigma_beta, sigma_a
         ), call. = FALSE)
  }
  for (sweep_index in seq_len(warmup)) {
    warm_state <- tryCatch(
      rat_dpglm_state(
        beta, random_effects, weights, support, X, exposure, rat_index
      ),
      error = function(error) NULL
    )
    if (is.null(warm_state) || any(!is.finite(warm_state$theta)) ||
        any(!is.finite(warm_state$btheta))) break
    warm_log_u <- log(stats::rexp(n)) - warm_state$btheta
    warm_log_psi <- rat_dpglm_log_psi(warm_log_u, warm_state$theta, support)
    if (any(!is.finite(warm_log_u)) || any(!is.finite(warm_log_psi))) break
    weights <- rat_dpglm_sample_weights(shape, warm_log_psi, weight_floor)
    if (any(!is.finite(weights))) break
    warm_lp <- rat_dpglm_logposterior(
      beta, random_effects, weights, support, y, X, exposure, rat_index,
      mu_beta, sigma_beta, sigma_a
    )
    if (is.finite(warm_lp) && warm_lp > best_lp) {
      best_lp <- warm_lp
      best_weights <- weights
    }
  }
  weights <- best_weights
  current_lp <- best_lp
  # Exact log-weights for the conditional alpha update (initial values that sit
  # at the floor are replaced as soon as their block is updated).
  track_log_w <- !is.null(alpha_prior) && alpha_update == "conditional"
  log_w <- log(weights)

  state <- rat_dpglm_state(
    beta, random_effects, weights, support, X, exposure, rat_index
  )
  
  # u_i | w, theta ~ Exponential(rate = exp(b_i)); stored on the log scale.
  log_u <- log(stats::rexp(n)) - state$btheta
  
  beta_samples <- matrix(NA_real_, iter, p)
  a_samples <- matrix(NA_real_, iter, n_rats)
  theta_samples <- matrix(NA_real_, iter, n)
  u_samples <- matrix(NA_real_, iter, n)
  mu_o_samples <- vector("list", iter)
  Jstar_samples <- vector("list", iter)
  weight_samples <- matrix(NA_real_, iter, length(support))
  loglik_samples <- numeric(iter)
  sigma_a_samples <- rep(sigma_a, iter)
  alpha_samples <- rep(alpha, iter)
  tau_samples <- matrix(NA_real_, iter, length(beta_groups))
  alpha_accept <- 0L
  alpha_accept_post <- 0L
  g0_prob_pos <- pmax(g0_prob, .Machine$double.xmin)
  crm_log_r_samples <- rep(NA_real_, iter)
  crm_data_samples <- rep(NA_real_, iter)
  crm_mass_samples <- rep(NA_real_, iter)
  crm_laplace_samples <- rep(NA_real_, iter)
  beta_accept <- 0L
  a_accept <- 0L
  joint_accept <- 0L
  crm_accept <- 0L
  beta_accept_post <- 0L
  a_accept_post <- 0L
  joint_accept_post <- 0L
  crm_accept_post <- 0L
  post_iterations <- 0L
  rat_rows <- split(seq_len(n), factor(rat_index, levels = seq_len(n_rats)))
  # Design for xi = (beta, a): eta_i = X_i beta + a_r(i) (+ log exposure).
  Z_mala <- cbind(X, outer(rat_index, seq_len(n_rats), "==") * 1)
  residual_index <- which(nstar_support == 0L)
  store_state <- function(iteration) {
    beta_samples[iteration, ] <<- beta
    a_samples[iteration, ] <<- random_effects
    theta_samples[iteration, ] <<- state$theta
    u_samples[iteration, ] <<- exp(log_u)
    mu_o_samples[[iteration]] <<- list(
      locations = support[residual_index],
      jumps = weights[residual_index]
    )
    Jstar_samples[[iteration]] <<- weights[zstar + 1L]
    weight_samples[iteration, ] <<- weights
    loglik_samples[iteration] <<- rat_dpglm_loglik(y, state, weights)
    sigma_a_samples[iteration] <<- sigma_a
    alpha_samples[iteration] <<- alpha
    if (length(beta_groups)) {
      tau_samples[iteration, ] <<- vapply(beta_groups, function(g) sigma_beta[g[1]], numeric(1))
    }
  }
  store_state(1L)
  
  for (iteration in 2:iter) {
    beta_accepted_now <- 0L
    a_accepted_now <- 0L
    joint_accepted_now <- 0L

    if (MCMC_update == "fisher_mala") {
      # Joint (beta, a) Fisher-preconditioned MALA:
      #   xi* ~ N(xi + (eps^2/2) H^-1 g, eps^2 H^-1), H = Fisher info + prior
      # precision at xi; the reverse move uses g and H at xi*.
      current_joint <- c(beta, random_effects)
      cur_terms <- rat_dpglm_mala_terms(beta, random_effects, state, weights,
                                        support, y, Z_mala, mu_beta, sigma_beta, sigma_a)
      if (!is.null(cur_terms$R)) {
        drift <- backsolve(cur_terms$R, forwardsolve(t(cur_terms$R), cur_terms$grad))
        proposed_joint <- current_joint + 0.5 * mala_step^2 * drift +
          mala_step * backsolve(cur_terms$R, stats::rnorm(joint_dim))
        beta_proposal <- proposed_joint[seq_len(p)]
        a_proposal <- proposed_joint[p + seq_len(n_rats)]
        prop_state <- tryCatch(
          rat_dpglm_state(beta_proposal, a_proposal, weights, support, X,
                          exposure, rat_index,
                          if (theta_warm_start) state$theta else NULL),
          error = function(error) NULL)
        proposal_lp <- if (is.null(prop_state) || any(!is.finite(prop_state$theta)) ||
                           any(!is.finite(prop_state$btheta))) -Inf else
          rat_dpglm_loglik(y, prop_state, weights) +
          sum(dnorm(beta_proposal, mu_beta, sigma_beta, log = TRUE)) +
          sum(dnorm(a_proposal, 0, sigma_a, log = TRUE))
        if (is.finite(proposal_lp)) {
          prop_terms <- rat_dpglm_mala_terms(beta_proposal, a_proposal, prop_state,
                                             weights, support, y, Z_mala, mu_beta,
                                             sigma_beta, sigma_a)
          if (!is.null(prop_terms$R)) {
            log_r_mala <- proposal_lp - current_lp +
              rat_dpglm_mala_logq(current_joint, proposed_joint, prop_terms, mala_step) -
              rat_dpglm_mala_logq(proposed_joint, current_joint, cur_terms, mala_step)
            if (is.finite(log_r_mala) && log(stats::runif(1)) < log_r_mala) {
              beta <- beta_proposal
              random_effects <- a_proposal
              current_lp <- proposal_lp
              joint_accept <- joint_accept + 1L
              joint_accepted_now <- 1L
              beta_accept <- beta_accept + 1L
              a_accept <- a_accept + n_rats
              beta_accepted_now <- 1L
              a_accepted_now <- n_rats
            }
          }
        }
      }

    } else if (!is.null(joint_chol)) {
      # Joint  (beta, a)  random-walk proposal:
      # xi* = xi + joint_step z' L, with L'L = joint_cov.
      # The proposal is symmetric, so proposal-density terms cancel.
      current_joint <- c(beta, random_effects)
      proposed_joint <- current_joint +
        joint_step * drop(stats::rnorm(joint_dim) %*% joint_chol)

      beta_proposal <- proposed_joint[seq_len(p)]
      a_proposal <- proposed_joint[p + seq_len(n_rats)]

      proposal_lp <- rat_dpglm_logposterior(
        beta_proposal, a_proposal, weights, support, y, X, exposure,
        rat_index, mu_beta, sigma_beta, sigma_a
      )

      if (is.finite(proposal_lp) &&
          log(stats::runif(1)) < proposal_lp - current_lp) {
        beta <- beta_proposal
        random_effects <- a_proposal
        current_lp <- proposal_lp

        joint_accept <- joint_accept + 1L
        joint_accepted_now <- 1L

        # Acceptance counts for beta and a.
        beta_accept <- beta_accept + 1L
        a_accept <- a_accept + n_rats
        beta_accepted_now <- 1L
        a_accepted_now <- n_rats
      }

    } else {
      # Random-walk update: beta as one block, then each a_r separately.
      beta_proposal <- if (is.null(beta_chol)) {
        beta + stats::rnorm(p, 0, beta_step)
      } else {
        beta + beta_step * drop(stats::rnorm(p) %*% beta_chol)
      }

      proposal_lp <- rat_dpglm_logposterior(
        beta_proposal, random_effects, weights, support, y, X, exposure,
        rat_index, mu_beta, sigma_beta, sigma_a
      )

      if (is.finite(proposal_lp) &&
          log(stats::runif(1)) < proposal_lp - current_lp) {
        beta <- beta_proposal
        current_lp <- proposal_lp
        beta_accept <- beta_accept + 1L
        beta_accepted_now <- 1L
      }

      for (rat in seq_len(n_rats)) {
        rows <- rat_rows[[rat]]
        if (!length(rows)) next

        current_block <- rat_dpglm_block_loglik(
          random_effects[rat], beta, rows, weights, support, X, exposure, y
        ) + dnorm(random_effects[rat], 0, sigma_a, log = TRUE)

        proposed_value <- random_effects[rat] + stats::rnorm(1, 0, a_step)

        proposed_block <- rat_dpglm_block_loglik(
          proposed_value, beta, rows, weights, support, X, exposure, y
        ) + dnorm(proposed_value, 0, sigma_a, log = TRUE)

        if (is.finite(proposed_block) && is.finite(current_block) &&
            log(stats::runif(1)) < proposed_block - current_block) {
          random_effects[rat] <- proposed_value
          a_accept <- a_accept + 1L
          a_accepted_now <- a_accepted_now + 1L
        }
      }
    }

    # Robbins-Monro tuning. Adaptation stops at adapt_until.  The joint
    # multivariate proposal targets acceptance 0.234.  Without joint_cov, retain
    # the beta-block and scalar-random-effect adaptation.
    if (iteration <= adapt_until) {
      gain <- 1 / sqrt(iteration)
      if (MCMC_update == "fisher_mala") {
        mala_step <- min(max(
          mala_step * exp(gain * (joint_accepted_now - 0.574)), 1e-3), 3)
      } else if (!is.null(joint_chol)) {
        joint_step <- min(max(
          joint_step * exp(gain * (joint_accepted_now - 0.234)), 1e-6), 1e3)
      } else {
        beta_step <- min(max(
          beta_step * exp(gain * (beta_accepted_now - 0.234)), 1e-6), 1e3)
        a_step <- min(max(
          a_step * exp(gain * (a_accepted_now / n_rats - 0.44)), 1e-6), 1e3)
      }
    }


    # Conjugate (Gibbs) updates of the variance hyperparameters, if requested.
    # Half-Cauchy(0, A) on a scale s:  s^2 | xi ~ IG(1/2, 1/xi), xi ~ IG(1/2, 1/A^2),
    # so  s^2 | b, xi ~ IG((m + 1)/2, 1/xi + sum(b^2)/2) for m coefficients b with
    # N(0, s^2) priors, and  xi | s^2 ~ IG(1, 1/A^2 + 1/s^2).
    if (sigma_a_prior == "halfcauchy") {
      sigma_a <- sqrt(1 / stats::rgamma(1, (n_rats + 1) / 2,
                                         rate = 1 / xi_sigma_a + sum(random_effects^2) / 2))
      xi_sigma_a <- 1 / stats::rgamma(1, 1, rate = 1 / sigma_a_scale^2 + 1 / sigma_a^2)
    } else if (sigma_a_prior == "invgamma") {
      sigma_a <- sqrt(1 / stats::rgamma(1, sigma_a_ig[1] + n_rats / 2,
                                         rate = sigma_a_ig[2] + sum(random_effects^2) / 2))
    }
    for (gi in seq_along(beta_groups)) {
      g <- beta_groups[[gi]]
      tau2 <- 1 / stats::rgamma(1, (length(g) + 1) / 2,
                                rate = 1 / xi_beta_groups[gi] + sum(beta[g]^2) / 2)
      xi_beta_groups[gi] <- 1 / stats::rgamma(1, 1, rate = 1 / beta_group_scale^2 + 1 / tau2)
      sigma_beta[g] <- sqrt(tau2)
    }

    state <- rat_dpglm_state(
      beta, random_effects, weights, support, X, exposure, rat_index,
      if (theta_warm_start) state$theta else NULL
    )
    
    # Exact draw of the auxiliary variables given the current weights.
    log_u <- log(stats::rexp(n)) - state$btheta

    # Optional exact rescaling move for the total CRM mass (see crm_mass_move).
    # Under (w, u) -> (c w, u / c) every u_i T_i(w), theta and the normalized
    # weights are invariant; the remaining density of M = sum(w) is
    # M^(alpha - 1) exp(-M), so M* ~ Gamma(alpha, 1) is a Gibbs draw.
    if (isTRUE(crm_mass_move)) {
      log_c <- log(stats::rgamma(1, shape = alpha, rate = 1)) - log(sum(weights))
      weights <- pmax(weights * exp(log_c), weight_floor)
      log_w <- log_w + log_c
      log_u <- log_u - log_c
      state$btheta <- state$btheta + log_c
    }
    
    crm_accepted_now <- 0L
    if (crm_blocks > 1L) {
      # Block CRM update: random partition of the support; each block gets the
      # conjugate proposal given (u, theta) and an exact MH correction, with u
      # held fixed (a valid update of the u-augmented target).
      blocks <- split(sample.int(length(support)),
                      rep_len(seq_len(crm_blocks), length(support)))
      for (block in blocks) {
        log_psi <- rat_dpglm_log_psi(log_u, state$theta, support)
        proposed_weights <- weights
        if (track_log_w) {
          proposed_log_w <- log_w
          proposed_log_w[block] <- rat_dpglm_sample_log_weights(shape[block], log_psi[block])
          proposed_weights[block] <- pmax(exp(proposed_log_w[block]), weight_floor)
        } else {
          proposed_weights[block] <- rat_dpglm_sample_weights(
            shape[block], log_psi[block], weight_floor)
        }
        proposed_state <- tryCatch(
          rat_dpglm_state(beta, random_effects, proposed_weights, support, X,
                          exposure, rat_index,
                          if (theta_warm_start) state$theta else NULL),
          error = function(error) NULL)
        if (is.null(proposed_state) || any(!is.finite(proposed_state$theta)) ||
            any(!is.finite(proposed_state$btheta))) next
        log_psi_star <- rat_dpglm_log_psi(log_u, proposed_state$theta, support)
        log_r <- rat_dpglm_crm_block_log_ratio(
          y, log_u, shape, state$theta, proposed_state$theta, support, weights,
          proposed_weights, log_psi, log_psi_star, block)
        crm_log_r_samples[iteration] <- log_r
        crm_data_samples[iteration] <- attr(log_r, "data_term")
        crm_mass_samples[iteration] <- attr(log_r, "mass_term")
        crm_laplace_samples[iteration] <- attr(log_r, "laplace_term")
        if (is.finite(log_r) && log(stats::runif(1)) < log_r) {
          weights <- proposed_weights
          if (track_log_w) log_w <- proposed_log_w
          state <- proposed_state
          crm_accepted_now <- crm_accepted_now + 1 / crm_blocks
        }
      }
      crm_accept <- crm_accept + crm_accepted_now
    } else {
    log_psi <- rat_dpglm_log_psi(log_u, state$theta, support)
    if (track_log_w) {
      proposed_log_w <- rat_dpglm_sample_log_weights(shape, log_psi)
      proposed_weights <- pmax(exp(proposed_log_w), weight_floor)
    } else {
      proposed_weights <- rat_dpglm_sample_weights(shape, log_psi, weight_floor)
    }
    proposed_state <- tryCatch(
      rat_dpglm_state(
        beta, random_effects, proposed_weights, support, X, exposure, rat_index,
        if (theta_warm_start) state$theta else NULL
      ),
      error = function(error) NULL
    )
    if (!is.null(proposed_state) && all(is.finite(proposed_state$theta)) &&
        all(is.finite(proposed_state$btheta))) {
      log_psi_star <- rat_dpglm_log_psi(
        log_u, proposed_state$theta, support
      )
      log_r <- rat_dpglm_crm_log_ratio(
        y = y,
        log_u = log_u,
        shape = shape,
        theta = state$theta,
        theta_star = proposed_state$theta,
        support = support,
        weights = weights,
        proposed_weights = proposed_weights,
        log_psi = log_psi,
        log_psi_star = log_psi_star,
        btheta_current = state$btheta,
        btheta_star_proposed = proposed_state$btheta
      )
      crm_log_r_samples[iteration] <- log_r
      crm_data_samples[iteration] <- attr(log_r, "data_term")
      crm_mass_samples[iteration] <- attr(log_r, "mass_term")
      crm_laplace_samples[iteration] <- attr(log_r, "laplace_term")
      if (is.finite(log_r) && log(stats::runif(1)) < log_r) {
        weights <- proposed_weights
        if (track_log_w) log_w <- proposed_log_w
        state <- proposed_state
        crm_accept <- crm_accept + 1L
        crm_accepted_now <- 1L
      }
    }
    }

    # Joint update of alpha and all CRM weights (only if alpha_prior is given).
    # Propose log alpha* = log alpha + alpha_step * z, then w*_k ~
    # Gamma(alpha* G0(k) + n_k, 1 + psi(k)) given the current (u, theta).  The
    # log w terms cancel exactly, leaving
    #   log R = log p(alpha*) - log p(alpha) + log(alpha*/alpha)
    #         + sum_i (theta*_i - theta_i) y_i
    #         + sum_i u_i [T_i(theta; w* + w) - T_i(theta*; w* + w)]
    #         + sum_k [c_k log(1 + psi*(k)) - c*_k log(1 + psi(k))]
    #         + sum_k [lgamma(c*_k) - lgamma(c_k)
    #                  - lgamma(alpha* G0(k)) + lgamma(alpha G0(k))],
    # which reduces to the CRM ratio when alpha* = alpha.
    alpha_accepted_now <- 0L
    if (track_log_w) {
      # Slice sampler (stepping out, then shrinkage) for z = log(alpha) from
      #   log p(z | w) = log p(alpha) + z + alpha * sum_k G0(k) log w_k
      #                  - sum_k lgamma(alpha G0(k)) + const.
      sum_g_logw <- sum(g0_prob_pos * log_w)
      log_target <- function(z) {
        a <- exp(z)
        stats::dgamma(a, alpha_prior[1], alpha_prior[2], log = TRUE) + z +
          a * sum_g_logw - sum(lgamma(a * g0_prob_pos))
      }
      z0 <- log(alpha); level <- log_target(z0) - stats::rexp(1)
      left <- z0 - stats::runif(1); right <- left + 1
      for (step in 1:50) { if (log_target(left) < level) break; left <- left - 1 }
      for (step in 1:50) { if (log_target(right) < level) break; right <- right + 1 }
      repeat {
        z1 <- stats::runif(1, left, right)
        if (log_target(z1) >= level) break
        if (z1 < z0) left <- z1 else right <- z1
      }
      alpha <- exp(z1)
      shape <- alpha * g0_prob + nstar_support
      alpha_accept <- alpha_accept + 1L
      alpha_accepted_now <- 1L
    }
    if (!is.null(alpha_prior) && alpha_update == "joint") {
      log_psi <- rat_dpglm_log_psi(log_u, state$theta, support)
      alpha_star <- alpha * exp(alpha_step * stats::rnorm(1))
      shape_star <- alpha_star * g0_prob + nstar_support
      proposed_weights <- rat_dpglm_sample_weights(shape_star, log_psi, weight_floor)
      proposed_state <- tryCatch(
        rat_dpglm_state(beta, random_effects, proposed_weights, support, X,
                        exposure, rat_index,
                        if (theta_warm_start) state$theta else NULL),
        error = function(error) NULL)
      if (!is.null(proposed_state) && all(is.finite(proposed_state$theta)) &&
          all(is.finite(proposed_state$btheta))) {
        log_psi_star <- rat_dpglm_log_psi(log_u, proposed_state$theta, support)
        mixed <- proposed_weights + weights
        log_r_alpha <-
          stats::dgamma(alpha_star, alpha_prior[1], alpha_prior[2], log = TRUE) -
          stats::dgamma(alpha, alpha_prior[1], alpha_prior[2], log = TRUE) +
          log(alpha_star) - log(alpha) +
          sum((proposed_state$theta - state$theta) * y) +
          sum(exp(log_u + rat_dpglm_btheta(state$theta, support, mixed)) -
                exp(log_u + rat_dpglm_btheta(proposed_state$theta, support, mixed))) +
          sum(shape * rat_dpglm_log1p_exp(log_psi_star) -
                shape_star * rat_dpglm_log1p_exp(log_psi)) +
          sum(lgamma(shape_star) - lgamma(shape) -
                lgamma(alpha_star * g0_prob_pos) + lgamma(alpha * g0_prob_pos))
        if (is.finite(log_r_alpha) && log(stats::runif(1)) < log_r_alpha) {
          alpha <- alpha_star
          shape <- shape_star
          weights <- proposed_weights
          state <- proposed_state
          alpha_accept <- alpha_accept + 1L
          alpha_accepted_now <- 1L
        }
      }
      if (iteration <= adapt_until) {
        alpha_step <- min(max(
          alpha_step * exp((alpha_accepted_now - 0.25) / sqrt(iteration)), 1e-3), 5)
      }
    }

    current_lp <- rat_dpglm_logposterior(
      beta, random_effects, weights, support, y, X, exposure, rat_index,
      mu_beta, sigma_beta, sigma_a,
      if (theta_warm_start) state$theta else NULL
    )
    
    if (iteration > adapt_until) {
      post_iterations <- post_iterations + 1L
      beta_accept_post <- beta_accept_post + beta_accepted_now
      a_accept_post <- a_accept_post + a_accepted_now
      joint_accept_post <- joint_accept_post + joint_accepted_now
      crm_accept_post <- crm_accept_post + crm_accepted_now
      alpha_accept_post <- alpha_accept_post + alpha_accepted_now
    }

    store_state(iteration)
    if (isTRUE(report_progress) && iteration %% max(1L, floor(iter / 10L)) == 0L) {
      message(sprintf("rat DPGLM iteration %d/%d", iteration, iter))
    }
  }
  
  list(
    model = "fcDPGLM",
    beta = beta_samples,
    a = a_samples,
    theta = theta_samples,
    u = u_samples,
    zstar = zstar,
    nstar = nstar,
    mu_o = mu_o_samples,
    Jstar = Jstar_samples,
    weights = weight_samples,
    support = support,
    shape = shape,
    g0 = list(family = g0, mu = g0_mu, size = g0_size, probability = g0_prob),
    sigma_a = sigma_a,
    sigma_a_samples = sigma_a_samples,
    alpha_samples = alpha_samples,
    tau_samples = tau_samples,
    loglik = loglik_samples,
    crm_log_r = crm_log_r_samples,
    crm_data = crm_data_samples,
    crm_mass = crm_mass_samples,
    crm_laplace = crm_laplace_samples,
    tuning = list(
      beta_step = beta_step,
      a_step = a_step,
      joint_step = joint_step,
      beta_cov_used = !is.null(beta_chol),
      joint_cov_used = !is.null(joint_chol),
      joint_dim = if (!is.null(joint_chol)) joint_dim else NULL,
      adapt_until = adapt_until,
      warmup = warmup,
      MCMC_update = MCMC_update,
      mala_step = if (MCMC_update == "fisher_mala") mala_step else NULL,
      crm_blocks = crm_blocks,
      sigma_a_prior = sigma_a_prior,
      alpha_prior = alpha_prior,
      alpha_step = alpha_step,
      beta_groups = beta_groups,
      theta_warm_start = theta_warm_start,
      use_cpp_kernels = isTRUE(use_cpp_kernels),
      crm_mass_move = isTRUE(crm_mass_move),
      alpha_update = alpha_update
    ),
    acceptance = list(
      beta = beta_accept / iter,
      random_effects = a_accept / (iter * n_rats),
      joint = if (!is.null(joint_chol) || MCMC_update == "fisher_mala")
        joint_accept / iter else NA_real_,
      crm = crm_accept / iter,
      post_adaptation = list(
        iterations = post_iterations,
        beta = if (post_iterations > 0L) beta_accept_post / post_iterations else NA_real_,
        random_effects = if (post_iterations > 0L)
          a_accept_post / (post_iterations * n_rats) else NA_real_,
        joint = if ((!is.null(joint_chol) || MCMC_update == "fisher_mala") &&
                    post_iterations > 0L)
          joint_accept_post / post_iterations else NA_real_,
        crm = if (post_iterations > 0L) crm_accept_post / post_iterations else NA_real_,
        alpha = if (!is.null(alpha_prior) && post_iterations > 0L)
          alpha_accept_post / post_iterations else NA_real_
      )
    )
  )
}