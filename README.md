# Nutrimumouse sCCA Train–Test Tutorial

A reproducible teaching workflow for sparse canonical correlation analysis (sCCA) using the `Nutrimumouse` dataset from `mixOmics`.

The analysis links:

- **X block:** 120 gene-expression variables
- **Y block:** 21 lipid variables
- **Samples:** 40 mice

## Workflow

### 1. Fit the sCCA model

` sCCA_mouse_3comp_train_test(6).R `

- Splits the data into 75% training and 25% held-out test sets
- Estimates scaling parameters and canonical weights using training data only
- Fits three sequential sCCA components
- Uses projection deflation for Components 2 and 3
- Applies training weights to the test data
- Exports weights, canonical scores, correlations, loadings, cross-loadings, selected variables, and sample metadata

### 2. Build interpretation outputs

`01_build_scca_interpretation_outputs.R`

- Reads existing sCCA outputs without refitting the model
- Performs automated quality-control checks
- Combines weights, loadings, cross-loadings, and selection summaries
- Produces downstream interpretation and validation tables

### 3. Generate the Figure 2-inspired visualization

`02_plot_fig2_inspired_mouse_scca.R`

- Uses the fixed train/test model outputs
- Visualizes train and test canonical correlations and scores
- Displays selected gene weights
- Reconstructs gene EV and FEV from training X-to-Y-score cross-loadings
- Displays selected lipid profiles using absolute Y-side weights
- Exports publication-ready PDF, SVG, and PNG figures

## Requirements

- R
- `mixOmics`
- `nscancor`
- `ggplot2`
- `dplyr`
- `tidyr`
- `patchwork`
- `scales`
- `ragg`
- `svglite`

## Important note

This repository is a **methods-inspired teaching reconstruction**, not an exact reproduction of the full Xu et al. (2023) analysis pipeline.

The current implementation uses fixed `keepX` values to control sparsity and applies paper-like test deflation for later components. The EV/FEV panels are reconstructed from squared cross-loadings and should not be interpreted as an exact implementation of the original paper’s non-sCCA FEV analysis.
