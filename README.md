# fcDPGLM: code and synthetic data

Code and a synthetic data set for the paper

> Alam, E., Juarez-Colunga, E., González, M. I. and Müller, P.
> *A Semiparametric Bayesian Generalized Linear Mixed Model with a Nonparametric
> Count Distribution for Seizure Data.* Submitted to *Biometrics*.

The repository contains the fcDPGLM sampler, the scripts for the seizure data
analysis in Section 5 of the paper, and their results on the synthetic data.

## Synthetic data

The seizure data analyzed in the paper are not publicly available. To allow the
analysis to be run, `data/synthetic_seizure_counts.csv` contains
simulated daily seizure counts with the same structure: 21 rats, two stages
and 618 rat-days.

| Column     | Description                                                          |
|------------|----------------------------------------------------------------------|
| `rat`      | rat identifier, 1-21                                                 |
| `stage`    | 1 = early period, 2 = chronic period                                 |
| `day`      | day within the stage (1, 2, ...)                                     |
| `exposure` | observed recording time on that day, in days (between 0 and 1)       |
| `y`        | number of seizures on that day                                       |

The days and exposures follow the observation schedule of the original study.
The counts were generated from the semiparametric model of the paper, with
mean curves, rat effects and reference distribution estimated from the
seizure data. This is the semiparametric setting with shape-difference mean
curves and 21 rats in Section 4. The rat effects were resampled with
replacement, so rat `r` in the synthetic data does not correspond to rat `r`
in the original study.
Summary statistics are close to those of the real data (mean 11.7 against
10.2 seizures per day, 25% against 29% zero counts), but results obtained
from the synthetic data are not the results reported in the paper.

## Repository layout

```
data/       synthetic_seizure_counts.csv
R/          fcdpglm.R: the fcDPGLM sampler (fit_rat_dpglm())
src/        C++ code for the tilting-parameter solve and likelihood kernels
scripts/    analysis scripts, run in order
results/    output of run_all.sh on the synthetic data
run_all.sh  runs all scripts
```

| Script                         | Output                                                           |
|--------------------------------|------------------------------------------------------------------|
| `01_fit_fcdpglm.R`             | fcDPGLM fit, four chains (`results/fcdpglm_fit.rds`, not tracked; see below) |
| `02_fit_glmms.R`               | Poisson and negative-binomial GLMMs and the model comparison (`results/glmm_comparison/`) |
| `03_posterior_summaries.R`     | posterior summaries, stage contrasts, rat effects, predictive checks (`results/fcdpglm/`) |
| `04_figures.R`                 | mean functionals and model comparison figures (`results/figures/`) |
| `05_count_distribution.R`      | dispersion index, zero and exceedance probabilities as functions of the mean (`results/figures/`) |
| `06_trace_plots.R`             | trace plots by chain (`results/figures/`)                        |

## Running the analysis

From the repository root:

```bash
./run_all.sh
```

or run the scripts one at a time, for example

```bash
Rscript scripts/01_fit_fcdpglm.R chains=4 iter=30000 adapt=5000 thin=5
```

The fcDPGLM fit runs the four chains in parallel (via `parallel::mclapply`,
so on Windows they run one after another); each chain took about 7.5 minutes,
and `run_all.sh` about 8 minutes in total, on a laptop with four free cores.
A short test run takes under a minute:

```bash
Rscript scripts/01_fit_fcdpglm.R chains=2 iter=600 adapt=300 thin=5
```

The fitted object `results/fcdpglm_fit.rds` (about 200 MB) is not included in
the repository; rerun `01_fit_fcdpglm.R` to recreate it.

## Requirements

R (tested with 4.5.1) and the packages `Rcpp`, `RcppArmadillo`, `glmmTMB`,
`coda`, `MASS`, `ggplot2`, `patchwork`, `latex2exp`, `dplyr` and `tidyr`.
The C++ code is compiled on first use with `Rcpp::sourceCpp()`, which needs a
C++ compiler (Xcode command line tools on macOS, Rtools on Windows).
`results/sessionInfo.txt` lists the package versions used.

## License

The code is released under the MIT License (see `LICENSE`).

## Contact

Entejar Alam, entejar@utexas.edu
