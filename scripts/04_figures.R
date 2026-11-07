#!/usr/bin/env Rscript
# ============================================================================
# Figures for the seizure count analysis (run from the repository root):
#   fig_case_functionals.pdf  (A) stage-specific mean curves m_d(t),
#                             (B) Delta(t) with pointwise and simultaneous bands,
#                             (C) shape contrast S(t) with simultaneous band
#   fig_case_comparison.pdf   (A) m_d(t) and (B) Delta(t) under fcDPGLM and the
#                             Poisson and NB GLMMs, (C) predictive reproduction
#   rat_random_effects.pdf    exp(a_r) by rat
# Inputs: the fcDPGLM fit (scripts/01_fit_fcdpglm.R) and the GLMM comparison
#   (scripts/02_fit_glmms.R).
# Usage: Rscript scripts/04_figures.R [fit.rds] [comparison_dir] [out_dir]
# ============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(splines); library(patchwork) })
TeX <- latex2exp::TeX
args <- commandArgs(trailingOnly = TRUE)
fit_path <- if (length(args) >= 1) args[1] else "results/fcdpglm_fit.rds"
cmp_dir  <- if (length(args) >= 2) args[2] else "results/glmm_comparison"
out_dir  <- if (length(args) >= 3) args[3] else "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

o <- readRDS(fit_path); f <- o$rat_dpglm; dat <- as.data.frame(o$data); X <- o$design_matrix
cn <- colnames(X); b <- f$beta
basis <- ns(dat$t, df = 3, intercept = FALSE)
tod <- function(day, s) { d <- dat$day[dat$stage == s]; (day - mean(d)) / sd(d) }
rng <- sapply(split(dat$day, dat$stage), range); grid <- seq(max(rng[1, ]), min(rng[2, ]))
B1 <- predict(basis, tod(grid, 1)); B2 <- predict(basis, tod(grid, 2))
i1 <- grep("^stage1_spline", cn); i2 <- grep("^stage2_spline", cn)
e1 <- b[, i1] %*% t(B1) + b[, "intercept"]
e2 <- b[, i2] %*% t(B2) + b[, "intercept"] + b[, "stage2"]
m1 <- exp(e1); m2 <- exp(e2); dl <- m2 - m1
sh <- (e2 - rowMeans(e2)) - (e1 - rowMeans(e1))

q <- function(M, p) apply(M, 2, quantile, p)
sim_band <- function(M) {                     # 95% simultaneous (sup-t) band
  mu <- colMeans(M); s <- apply(M, 2, sd)
  crit <- quantile(apply(abs(sweep(M, 2, mu)) / rep(s, each = nrow(M)), 1, max), .95)
  list(lo = mu - crit * s, hi = mu + crit * s)
}
# Figures are saved at (about) their printed width so text prints at full size.
th <- theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                     legend.key.size = grid::unit(0.8, "lines"))
xlab_t <- TeX("Day within stage, $t$")

# ---- Figure: case-study functionals ----------------------------------------
cur <- rbind(data.frame(day = grid, stage = "Stage 1", mean = colMeans(m1), lo = q(m1, .025), hi = q(m1, .975)),
             data.frame(day = grid, stage = "Stage 2", mean = colMeans(m2), lo = q(m2, .025), hi = q(m2, .975)))
pA <- ggplot(cur, aes(day, mean, colour = stage, fill = stage)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, colour = NA) + geom_line(linewidth = .8) +
  scale_colour_manual(values = c("#1b6ca8", "#d1495b")) + scale_fill_manual(values = c("#1b6ca8", "#d1495b")) +
  labs(x = xlab_t, y = TeX("$m_d(t)$ (seizures per day)"), colour = NULL, fill = NULL) + th
sb <- sim_band(dl)
dd <- data.frame(day = grid, mean = colMeans(dl), lo = q(dl, .025), hi = q(dl, .975), slo = sb$lo, shi = sb$hi)
pB <- ggplot(dd, aes(day)) +
  geom_ribbon(aes(ymin = slo, ymax = shi, fill = "Simultaneous 95%"), alpha = .25) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "Pointwise 95%"), alpha = .35) +
  geom_line(aes(y = mean), linewidth = .8) + geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("Pointwise 95%" = "grey45", "Simultaneous 95%" = "grey75")) +
  labs(x = xlab_t, y = TeX("$\\Delta(t) = m_2(t) - m_1(t)$"), fill = NULL) + th
sbS <- sim_band(sh)
ds <- data.frame(day = grid, mean = colMeans(sh), slo = sbS$lo, shi = sbS$hi)
pC <- ggplot(ds, aes(day)) +
  geom_ribbon(aes(ymin = slo, ymax = shi, fill = "Simultaneous 95%"), alpha = .4) +
  geom_line(aes(y = mean), linewidth = .8) + geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("Simultaneous 95%" = "grey75")) +
  labs(x = xlab_t, y = TeX("Shape contrast $S(t)$"), fill = NULL) + th
ggsave(file.path(out_dir, "fig_case_functionals.pdf"),
       (pA | pB | pC) + plot_annotation(tag_levels = "A"), width = 7.5, height = 3.2)

# ---- Figure: comparison with parametric mixed models ------------------------
mlev <- c("fcDPGLM", "Poisson GLMM", "NB GLMM")
mlab <- parse(text = c("italic(fc)*DPGLM", "'Poisson GLMM'", "'NB GLMM'"))
mcol <- c("fcDPGLM" = "#1b6ca8", "Poisson GLMM" = "#8a8a8a", "NB GLMM" = "#d1495b")
sc <- read.csv(file.path(cmp_dir, "full_data_stage_curves_and_delta.csv"))
sc$model <- factor(sc$model, mlev)
cc <- rbind(data.frame(sc[, c("model", "day")], stage = "Stage 1", mean = sc$stage1_mean, lo = sc$stage1_lower, hi = sc$stage1_upper),
            data.frame(sc[, c("model", "day")], stage = "Stage 2", mean = sc$stage2_mean, lo = sc$stage2_lower, hi = sc$stage2_upper))
scale_m <- list(scale_colour_manual(values = mcol, labels = mlab), scale_fill_manual(values = mcol, labels = mlab),
                scale_linetype_manual(values = c(1, 2, 4), labels = mlab))
qA <- ggplot(cc, aes(day, mean, colour = model, fill = model, linetype = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .1, colour = NA) + geom_line(linewidth = .7) +
  facet_wrap(~ stage) + scale_m +
  labs(x = xlab_t, y = TeX("$m_d(t)$ (seizures per day)"), colour = NULL, fill = NULL, linetype = NULL) + th
qB <- ggplot(sc, aes(day, delta_mean, colour = model, fill = model, linetype = model)) +
  geom_ribbon(aes(ymin = delta_lower, ymax = delta_upper), alpha = .1, colour = NA) + geom_line(linewidth = .7) +
  geom_hline(yintercept = 0, linetype = 3) + scale_m +
  labs(x = xlab_t, y = TeX("$\\Delta(t)$"), colour = NULL, fill = NULL, linetype = NULL) + th
pp <- read.csv(file.path(cmp_dir, "full_data_predictive_checks.csv"))
pp <- pp[pp$statistic %in% c("zero_rate", "variance", "p90", "p95", "maximum"), ]
pp$model <- factor(pp$model, rev(mlev))
pp$stat <- factor(pp$statistic, rev(c("zero_rate", "variance", "p90", "p95", "maximum")),
                  rev(c("Zero proportion", "Variance", "90th percentile", "95th percentile", "Maximum")))
qC <- ggplot(pp, aes(predictive_median / observed, stat, colour = model)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = lower95 / observed, xmax = upper95 / observed),
                  position = position_dodge(width = .6), size = .25) +
  scale_colour_manual(values = mcol, labels = rev(mlab), breaks = rev(mlev)) +
  labs(x = "Predictive / observed", y = NULL, colour = NULL) + th + guides(colour = "none")
ggsave(file.path(out_dir, "fig_case_comparison.pdf"),
       ((qA / qB) | qC) + plot_layout(widths = c(1.5, 1), guides = "collect") +
         plot_annotation(tag_levels = "A") & theme(legend.position = "bottom"),
       width = 7.5, height = 5.6)

# ---- Rat-specific effects --------------------------------------------------
a <- f$a; rats <- unique(as.character(dat$int_id)); ea <- exp(a)
re <- data.frame(rat = rats, mean = colMeans(ea), lo = q(ea, .025), hi = q(ea, .975))
re$rat <- factor(re$rat, re$rat[order(re$mean)])
pr <- ggplot(re, aes(mean, rat)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y") + geom_point(size = 1.8) +
  labs(x = TeX("Multiplicative rat effect $\\exp(a_r)$"), y = "Rat") + th
ggsave(file.path(out_dir, "rat_random_effects.pdf"), pr, width = 5.5, height = 5.5)
cat("Wrote fig_case_functionals.pdf, fig_case_comparison.pdf, rat_random_effects.pdf to", out_dir, "\n")
