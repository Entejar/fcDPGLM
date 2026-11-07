#!/usr/bin/env Rscript
# Trace plots by chain for the fcDPGLM fit (run from the repository root):
#   Rscript scripts/06_trace_plots.R [fit.rds] [out_dir]
suppressPackageStartupMessages({ library(ggplot2); library(splines) })
args <- commandArgs(trailingOnly = TRUE)
fit_path <- if (length(args) >= 1) args[1] else "results/fcdpglm_fit.rds"
out_dir  <- if (length(args) >= 2) args[2] else "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
o <- readRDS(fit_path); f <- o$rat_dpglm; dat <- as.data.frame(o$data); X <- o$design_matrix; cn <- colnames(X)
ch <- factor(f$chain_id); it <- ave(seq_along(ch), ch, FUN = seq_along) * f$tuning$thin
# Mean functionals
basis <- ns(dat$t, df = 3, intercept = FALSE)
tod <- function(day, s) { d <- dat$day[dat$stage == s]; (day - mean(d)) / sd(d) }
rng <- sapply(split(dat$day, dat$stage), range); grid <- seq(max(rng[1, ]), min(rng[2, ]))
B1 <- predict(basis, tod(grid, 1)); B2 <- predict(basis, tod(grid, 2)); Ba <- predict(basis, dat$t)
b <- f$beta; i1 <- grep("^stage1_spline", cn); i2 <- grep("^stage2_spline", cn)
e1 <- b[, i1] %*% t(B1) + b[, 1]; e2 <- b[, i2] %*% t(B2) + b[, 1] + b[, 2]
g1 <- b[, i1] %*% t(Ba) + b[, 1]; g2 <- b[, i2] %*% t(Ba) + b[, 1] + b[, 2]
long <- function(vals, labs) do.call(rbind, lapply(seq_along(labs), function(j)
  data.frame(iteration = it, chain = ch, value = vals[[j]],
             quantity = factor(labs[j], levels = labs))))
th <- theme_bw() + theme(legend.position = "bottom")
p1 <- ggplot(long(list(b[, 1], f$loglik, log(f$alpha_samples), f$sigma_a_samples, b[, 2]),
                  c("'Intercept'", "'Log likelihood'", "log~alpha", "sigma[a]", "'Stage effect'~beta[stage]")),
             aes(iteration, value, colour = chain)) + geom_line(linewidth = 0.15, alpha = 0.7) +
  facet_wrap(~ quantity, scales = "free_y", ncol = 1, labeller = label_parsed) +
  labs(x = "Iteration (after adaptation)", y = NULL, colour = "Chain") + th
ggsave(file.path(out_dir, "trace_parameters_by_chain.pdf"), p1, width = 7.5, height = 9)
p2 <- ggplot(long(list(rowMeans(exp(e2) - exp(e1)), rowMeans(e2 - e1), rowMeans(exp(g2)) - rowMeans(exp(g1))),
                  c("bar(Delta)", "'Average log rate ratio'", "'Rat-day average'~bar(Delta)")),
             aes(iteration, value, colour = chain)) + geom_line(linewidth = 0.15, alpha = 0.7) +
  facet_wrap(~ quantity, scales = "free_y", ncol = 1, labeller = label_parsed) +
  labs(x = "Iteration (after adaptation)", y = NULL, colour = "Chain") + th
ggsave(file.path(out_dir, "trace_functionals_by_chain.pdf"), p2, width = 7.5, height = 6)
cat("Wrote trace plots to", out_dir, "\n")
