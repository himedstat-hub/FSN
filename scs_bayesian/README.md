# Bayesian SCS prognosis analysis

This directory contains R/CmdStan code for the Bayesian hierarchical logistic-regression analyses reported in:

> Bayesian Hierarchical Modeling of Six-Month Prognosis After Permanent Spinal Cord Stimulation Implantation

## Scope

The model estimates six-month prognosis among patients who underwent permanent SCS implantation after a successful trial. It does not model trial success, progression to implantation, causal treatment effects, or longer-term outcomes.

No patient-level data, fitted model objects, posterior draws, or generated outputs are included in this repository.

## Software

The analysis was run with R 4.6.0, CmdStanR 0.9.0, and CmdStan 2.39.0. Required R packages are `cmdstanr` and `posterior`.

Install CmdStan and set its path before running the analysis. For example:

```r
cmdstanr::install_cmdstan(version = "2.39.0")
cmdstanr::set_cmdstan_path(cmdstanr::cmdstan_path())
```

## Input data

Provide a de-identified CSV file through the `SCS_DATA_CSV` environment variable, or place it at `data/analysis_cohort.csv`. The file must not be committed to GitHub.

Required variables are:

- `Implantation` (1 = permanent implantation)
- `6M_Success` (1 = at least 50% NRS reduction at 6 months)
- `Age`, `BMI`, `Duration`, `MED`
- `NRS` (or the legacy source-column name `PreVAS`, which contains the pre-implantation NRS score)
- `Sex`, `Alcohol`, `Smoking`, `Psy`
- `Pay`, `UpperLower`, `Dx`, and `Hospital`

For the final analytic cohort, the script checks for 246 permanently implanted patients, 121 six-month successes, and no missing candidate predictor values.

## Reproducing the primary analysis

```r
Sys.setenv(
  CMDSTAN = cmdstanr::cmdstan_path(),
  SCS_DATA_CSV = "path/to/deidentified_analysis_cohort.csv",
  SCS_R_OUT = "outputs"
)
Sys.setenv(SCS_STAGE = "primary")
source("analysis/fit_scs_bayesian.R")
```

The primary model uses diagnosis-specific random intercepts with partial pooling, weakly informative priors, four MCMC chains, 1,500 warm-up iterations, 1,500 post-warm-up iterations, and `adapt_delta = 0.99`.

## Additional stages

Set `SCS_STAGE` to one of the following:

- `sensitivity_all`
- `cv_all`, followed by `cv_performance`
- `cv_no_diagnosis_all`, followed by `cv_no_diagnosis_performance`
- `posterior_predictive_check`
- `prior_predictive_check`

The five-fold cross-validation stages refit the complete hierarchical model in each training fold. Scaling parameters are estimated from training data only.

## Repository structure

- `analysis/fit_scs_bayesian.R`: primary, sensitivity, posterior-predictive, and cross-validation analyses
- `stan/scs_hierarchical_logistic.stan`: primary model specification
- `.gitignore`: prevents data and generated outputs from being committed

