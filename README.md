# MUJHU

Code for analyses related to MUJHU projects.

## Repository structure
- `model/mrgsolve`: mrgsolve implementations of the published PK models used in this study

| File | Description |
|-----|-------------|
| `tdf_mapbayr.cpp` | TFV model used for MAP Bayesian estimation of individual parameters (EBEs) |
| `tdf.cpp` | TFV model used to simulate drug concentrations and derive exposure metrics |
| `3tc_mapbayr.cpp` | 3TC/FTC model used for MAP Bayesian estimation of EBEs |
| `3tc.cpp` | 3TC/FTC model used to simulate drug concentrations and derive exposure metrics |

References: Leung E., et al., 2023 CPT:PSP (PMID: 37814498)

- `script/`: analysis scripts
  - `metadata_da.R`: preprocessing of PK data and clinical covariates for TFV/3TC data from the Uganda study
  - `mb_da.R`: merge microbiome, metadata, and exposure metrics for TFV/3TC data from the Uganda study
  - `tfv/`: TFV analysis pipeline for the Uganda study
  - `3tc/`: 3TC analysis pipeline for the Uganda study (same workflow as TFV)
  - `mn/`: FTC analysis pipeline for the Minnesota study (same workflow as TFV)

Within each drug-specific folder (`tfv/`, `3tc/`, `mn/`), the following scripts are used:

| Script | Description |
|---|---|
| `mapbayes_XXX.R` | MAP Bayesian estimation of individual PK parameters using observed drug concentration|
| `ebesim_XXX.R` | Simulation of drug concentration–time profiles and derivation of exposure metrics|
| `eda_XXX.R` | Exploratory analysis of drug exposure metrics and microbiome features|
| `pcoa_XXX.R` | Overall microbiome diversity (alpha and beta diversity) impact evaluation across drug exposure groups|
| `heatmap_XXX.R` | Correlation heatmaps between microbial taxa and drug exposure|
| `diffa_XXX.R` | Differential abundance analysis between high and low exposure groups |
| `glmnet_XXX.R` | Regularized regression models (ridge and LASSO)|
| `randomforest_XXX.R` | Random forest model|
| `xgboost_XXX.R` | XGBoost model|
| `selection_XXX.R` | Feature selection across different machine learning methods|

The analysis workflow is identical across drugs, with the suffix XXX replaced by tfv, 3tc, or ftc_mn depending on the dataset.

- `renv.lock`: R package environment for reproducibility


## Reproducibility
This repository uses `renv` to manage the R environment.

```r
install.packages("renv")
renv::restore()
```

## Data availability
Individual-level data are not publicly available due to human subject privacy restrictions.  
Scripts are provided to reproduce analyses given access to the original datasets.

