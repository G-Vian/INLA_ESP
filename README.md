# Bayesian spatio-temporal model of dengue in São Paulo State (2010–2024)

Reproducible code for a Bayesian hierarchical spatio-temporal model of monthly
dengue incidence across the 645 municipalities of São Paulo State, Brazil,
fitted with integrated nested Laplace approximation (INLA). The model combines
distributed-lag non-linear terms (DLNM) for climate, a structured spatial
(Besag/ICAR) effect, global and region-specific seasonality, autoregressive
terms, and municipal socioeconomic indicators.

This repository reproduces the main model, the DLNM exposure–lag–response
surfaces, and the leakage-free rolling-origin predictive validation.

## Repository layout

```
.
├── R/
│   ├── main.R                          # orchestrator — run this
│   ├── 00_config.R                     # configuration, priors, paths, switches
│   ├── R00b_figures.R                  # publication figure theme
│   ├── 01_data.R                       # loads the processed cache, builds the panel
│   ├── 02_features_dlnm_lags.R         # cross-basis (DLNM) and autoregressive lags
│   ├── 03_model.R                      # stepwise model building and final INLA fit
│   ├── External_1_rolling_origin.R     # leakage-free rolling-origin validation
│   ├── External_2_predictive_intervals.R # out-of-sample predictive intervals + MAE
│   └── 05_diagnostics.R                # external-module dispatcher
├── data/                          # input data (see "Data" below — not versioned)
├── figures/                       # output figures (created at run time)
└── README.md
```

## Requirements

- R (≥ 4.2)
- [R-INLA](https://www.r-inla.org/download-install) (stable release)
- CRAN packages: `dplyr`, `tidyr`, `ggplot2`, `tibble`, `sf`, `dlnm`,
  `splines`, `pROC`, `patchwork`, `knitr`, `scales`, `here`

```r
install.packages(c("dplyr","tidyr","ggplot2","tibble","sf","dlnm","splines",
                   "pROC","patchwork","knitr","scales","here"))
install.packages("INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)
```

## Data

The pipeline reads a single processed cache file that bundles the analysis-ready
inputs (the municipal case panel, monthly climate, the socioeconomic indicators,
the municipal shapefile and the Regional Health Department lookup):

```
data/DADOS_BASE_MENSAL_V75.RData
```

Because of its size this file is distributed separately:

- **Download:** <ADD LINK — Google Drive / OneDrive / Zenodo>
- Place it in the `data/` directory before running the pipeline.

### Raw data sources

The raw data are public and were obtained from:

- **Dengue case counts** — Brazilian notifiable diseases system (SINAN),
  via the [Mosqlimate](https://mosqlimate.org/) data platform.
- **Sociodemographic indicators** — Brazilian Institute of Geography and
  Statistics (IBGE), 2010 Census.
- **Municipal boundaries** — IBGE municipal mesh (SP, 2024).

The scripts that turn the raw sources into the processed cache are not part of
this repository; the cache above is provided so that the model and its
validation can be reproduced directly.

## How to run

From the repository root:

```bash
Rscript R/main.R
```

The orchestrator runs the modules in order. On a first run it loads the cache,
builds the feature matrix, fits the model, runs the rolling-origin validation
and writes figures to `figures/`.

> **Note on run time.** The rolling-origin validation (`24_rolling_origin.R`)
> refits the model once per monthly origin (~140 fits). It is the slow step.
> It checkpoints each origin to `data/rolling_checkpoints/`, so an interrupted
> run resumes where it stopped. To reproduce only the main model and the DLNM
> surfaces, comment out the `24_` and `25_` lines in `R/main.R`.

## Validation design

Predictive performance is assessed with a **leakage-free rolling-origin**
scheme (module `24`): for each target month *t* the model is refitted using only
months up to *t − 1* and predicts month *t*. The autoregressive predictors for
month *t* are months *t − 1* and *t − 2*, both inside the training window, so no
information outside the training window enters any prediction. Module `25`
turns each origin's fitted values into genuine out-of-sample **posterior
predictive intervals** (Monte Carlo over the negative-binomial likelihood) and
reports the mean absolute error relative to the mean observed count.

## Model summary

- **Likelihood:** zero-inflated negative binomial (type 1).
- **Climate:** distributed-lag non-linear cross-bases for minimum temperature
  and mean precipitation.
- **Spatial:** Besag/ICAR effect on the municipal contiguity graph.
- **Temporal:** global cyclic seasonality, region-specific (DRS) cyclic
  seasonality, and two autoregressive terms on log-lagged cases.
- **Socioeconomic:** municipal deprivation and demographic indicators
  (2010 Census), entered per specification.
- **Offset:** log municipal population.
- **Priors:** penalised-complexity priors on the random-effect precisions.

## Citation

If you use this code, please cite the accompanying article:

> <ADD CITATION ONCE PUBLISHED>
