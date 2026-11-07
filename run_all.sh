#!/usr/bin/env bash
# Run the full analysis of data/synthetic_seizure_counts.csv (from the
# repository root).  The fcDPGLM fit runs four chains in parallel; with four
# free cores the whole script takes about 8 minutes.
set -euo pipefail
Rscript scripts/01_fit_fcdpglm.R chains=4 iter=30000 adapt=5000 thin=5
Rscript scripts/02_fit_glmms.R
Rscript scripts/03_posterior_summaries.R
Rscript scripts/04_figures.R
Rscript scripts/05_count_distribution.R
Rscript scripts/06_trace_plots.R
Rscript -e 'writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")'
