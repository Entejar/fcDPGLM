#!/usr/bin/env Rscript
# ============================================================================
# Conditional count distribution learned by fcDPGLM, as a function of the
# conditional mean mu (run from the repository root):
#   fig_case_dispersion.pdf  (A) dispersion index V(mu)/mu, (B) Pr(Y = 0 | mu),
#                            (C) exceedance probability Pr(Y >= 20 | mu)
# under fcDPGLM (posterior mean and 95% band, from the tilted reference
# measure of each retained draw) and the NB GLMM (V/mu = 1 + mu/size,
# Pr(Y = 0) = {size/(size + mu)}^size, band from the Wald interval for
# log(size); for the exceedance probability, quantiles over 2000 draws of
# log(size)); the Poisson GLMM gives V/mu = 1, exp(-mu) and its Poisson tail.
# The rug shows the posterior-mean conditional means m_i of the observed
# rat-days.
# Usage: Rscript scripts/05_count_distribution.R [fit.rds] [comparison_dir] [out_dir]
# ============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("R/fcdpglm.R")
TeX <- latex2exp::TeX
args <- commandArgs(trailingOnly = TRUE)
fit_path <- if (length(args) >= 1) args[1] else "results/fcdpglm_fit.rds"
cmp_dir  <- if (length(args) >= 2) args[2] else "results/glmm_comparison"
out_dir  <- if (length(args) >= 3) args[3] else "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rat_dpglm_ensure_getTheta()

o <- readRDS(fit_path); f <- o$rat_dpglm; dat <- as.data.frame(o$data)
k <- f$support
mu_grid <- exp(seq(log(0.25), log(40), length.out = 60))
draws <- round(seq(1, nrow(f$weights), length.out = 2000))

# Dispersion index, zero and exceedance probabilities of the tilted reference measure.
cut <- 20
DI <- P0 <- PC <- matrix(NA_real_, length(draws), length(mu_grid))
for (j in seq_along(draws)) {
  w <- f$weights[draws[j], ]; f0 <- w / sum(w)
  th <- getTheta_cpp(spt = k, f0 = f0, mu = mu_grid, thetaStart = NULL, logit = FALSE)$theta
  lp <- outer(th, k) + matrix(log(f0), length(mu_grid), length(k), byrow = TRUE)
  p <- exp(lp - apply(lp, 1, max)); p <- p / rowSums(p)
  m1 <- drop(p %*% k); v <- drop(p %*% k^2) - m1^2
  DI[j, ] <- v / m1; P0[j, ] <- p[, 1]; PC[j, ] <- rowSums(p[, k >= cut, drop = FALSE])
}
band <- function(M, model) data.frame(mu = mu_grid, model = model, est = colMeans(M),
                                      lo = apply(M, 2, quantile, .025), hi = apply(M, 2, quantile, .975))

nb <- readRDS(file.path(cmp_dir, "full_data_glmms.rds"))$negative_binomial
ls <- glmmTMB::fixef(nb)$disp[[1]]
ls_se <- sqrt(diag(vcov(nb, full = TRUE)))[["disp~(Intercept)"]]
size <- exp(ls + c(est = 0, lo = -1.96, hi = 1.96) * ls_se)
nb_di <- data.frame(mu = mu_grid, model = "NB GLMM", est = 1 + mu_grid / size[["est"]],
                    lo = 1 + mu_grid / size[["hi"]], hi = 1 + mu_grid / size[["lo"]])
p0nb <- function(s) (s / (s + mu_grid))^s
nb_p0 <- data.frame(mu = mu_grid, model = "NB GLMM", est = p0nb(size[["est"]]),
                    lo = pmin(p0nb(size[["lo"]]), p0nb(size[["hi"]])),
                    hi = pmax(p0nb(size[["lo"]]), p0nb(size[["hi"]])))
set.seed(1); sz <- exp(rnorm(2000, ls, ls_se))
PCnb <- t(sapply(sz, function(s) pnbinom(cut - 1, size = s, mu = mu_grid, lower.tail = FALSE)))
nb_pc <- data.frame(mu = mu_grid, model = "NB GLMM",
                    est = pnbinom(cut - 1, size = size[["est"]], mu = mu_grid, lower.tail = FALSE),
                    lo = apply(PCnb, 2, quantile, .025), hi = apply(PCnb, 2, quantile, .975))
po_pc <- data.frame(mu = mu_grid, model = "Poisson GLMM",
                    est = ppois(cut - 1, mu_grid, lower.tail = FALSE), lo = NA, hi = NA)
po_di <- data.frame(mu = mu_grid, model = "Poisson GLMM", est = 1, lo = NA, hi = NA)
po_p0 <- data.frame(mu = mu_grid, model = "Poisson GLMM", est = exp(-mu_grid), lo = NA, hi = NA)

# Posterior-mean conditional means of the observed rat-days.
rat_idx <- match(as.character(dat$int_id), unique(as.character(dat$int_id)))
m_obs <- dat$exposure * exp(drop(o$design_matrix %*% colMeans(f$beta)) + colMeans(f$a)[rat_idx])
rug <- data.frame(mu = m_obs[m_obs >= min(mu_grid) & m_obs <= max(mu_grid)])

mlev <- c("fcDPGLM", "Poisson GLMM", "NB GLMM")
mlab <- parse(text = c("italic(fc)*DPGLM", "'Poisson GLMM'", "'NB GLMM'"))
mcol <- c("fcDPGLM" = "#1b6ca8", "Poisson GLMM" = "#8a8a8a", "NB GLMM" = "#d1495b")
sc <- list(scale_colour_manual(values = mcol, labels = mlab, breaks = mlev),
           scale_fill_manual(values = mcol, labels = mlab, breaks = mlev),
           scale_linetype_manual(values = c("fcDPGLM" = 1, "Poisson GLMM" = 2, "NB GLMM" = 4),
                                 labels = mlab, breaks = mlev),
           scale_x_log10(breaks = c(0.25, 1, 2.5, 5, 10, 20, 40),
                         labels = c("0.25", "1", "2.5", "5", "10", "20", "40")))
th <- theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
panel <- function(d, ylab, log_y = FALSE) {
  g <- ggplot(d, aes(mu, est, colour = model, fill = model, linetype = model)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA, na.rm = TRUE) +
    geom_line(linewidth = .7) +
    geom_rug(data = rug, aes(x = mu), inherit.aes = FALSE, alpha = .25, length = grid::unit(.02, "npc")) +
    sc + labs(x = TeX("Conditional mean $\\mu$"), y = ylab, colour = NULL, fill = NULL, linetype = NULL) + th
  if (log_y) g + scale_y_log10() else g
}
dA <- rbind(band(DI, "fcDPGLM"), nb_di, po_di); dB <- rbind(band(P0, "fcDPGLM"), nb_p0, po_p0)
dC <- rbind(band(PC, "fcDPGLM"), nb_pc, po_pc)
dA$model <- factor(dA$model, mlev); dB$model <- factor(dB$model, mlev); dC$model <- factor(dC$model, mlev)
pA <- panel(dA, TeX("Dispersion index $V(\\mu)/\\mu$"), log_y = TRUE)
pB <- panel(dB, TeX("$\\Pr(Y=0\\,|\\,\\mu)$"))
pC <- panel(dC, TeX("$\\Pr(Y\\geq 20\\,|\\,\\mu)$"))
ggsave(file.path(out_dir, "fig_case_dispersion.pdf"),
       (pA | pB | pC) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
         theme(legend.position = "bottom"), width = 7.5, height = 3.0)

# Values at selected means.
at <- c(1, 2, 5, 10, 20, 30)
idx <- sapply(at, function(a) which.min(abs(mu_grid - a)))
tab <- data.frame(mu = round(mu_grid[idx], 2),
                  fc_DI = colMeans(DI)[idx], fc_DI_lo = apply(DI, 2, quantile, .025)[idx],
                  fc_DI_hi = apply(DI, 2, quantile, .975)[idx], nb_DI = nb_di$est[idx],
                  fc_P0 = colMeans(P0)[idx], nb_P0 = nb_p0$est[idx], pois_P0 = exp(-mu_grid[idx]),
                  fc_P20 = colMeans(PC)[idx], fc_P20_lo = apply(PC, 2, quantile, .025)[idx],
                  fc_P20_hi = apply(PC, 2, quantile, .975)[idx], nb_P20 = nb_pc$est[idx],
                  pois_P20 = po_pc$est[idx])
write.csv(tab, file.path(out_dir, "case_dispersion_summary.csv"), row.names = FALSE)
cat("NB size:", round(size, 3), "\n"); print(tab, digits = 3)
cat("Mode of fcDPGLM dispersion index at mu =", round(mu_grid[which.max(colMeans(DI))], 1), "\n")
cat("Wrote", file.path(out_dir, "fig_case_dispersion.pdf"), "\n")
