# Matrix autoregression estimators

The complete core MAR estimator code and documentation are in
[`docs/MAR_ESTIMATORS.md`](docs/MAR_ESTIMATORS.md):

- `fit_mar_pointmass()` estimates MAR(P) models with a point-mass spike-and-slab prior. Individual coefficient groups give CGS; supplied coefficient partitions give GGS.
- `fit_mar_ssvs()` estimates MAR(P) models with a continuous spike-and-slab prior and returns SSVS inclusion probabilities.

Both functions take a `G x N x T_train` array, estimate row and column coefficient matrices for each lag, and return sparse estimates, posterior inclusion probabilities (PIPs), covariance matrices, and optional MCMC draws. The document includes their Description, Usage, Inputs, and Output fields.

## R environment

Package requirements and installation instructions are in
[`docs/MAR_ESTIMATORS.md`](docs/MAR_ESTIMATORS.md#packages). To reproduce the
archived single-run outputs, start Windows R 4.6.1. Each experiment's data
generation block checks the required package versions and sets the reference
RNG and matrix-product options.

## Run the examples

In a fresh R session at the repository root, load the packages and run both
function definitions in `docs/MAR_ESTIMATORS.md`. In the relevant example
document, run these R code blocks in order:

1. Run **Data generation** to set the true parameters and generate `Y_train_list`.
2. Run **Core-function inputs** to fit the generated training data.

| Experiment | Example document | Methods |
| --- | --- | --- |
| Simulation 1: MAR(1) | [`docs/EXAMPLE_SIMULATION1.md`](docs/EXAMPLE_SIMULATION1.md) | CGS, SSVS |
| Simulation 2: MAR(3) | [`docs/EXAMPLE_SIMULATION2.md`](docs/EXAMPLE_SIMULATION2.md) | GGS, SSVS |
| Supplement S2: MAR(1) | [`docs/EXAMPLE_SUPPLEMENT_S2.md`](docs/EXAMPLE_SUPPLEMENT_S2.md) | CGS, GGS |

The first R block in each example document sets `n_reps <- 100L`. Change it
to `n_reps <- 1L` to run only the first replication. Fitted results are in
`fits_by_rep`.

## One-seed results

Each example document includes a full one-seed result section with
all `A_hat`, `B_hat`, `PIP_A`, `PIP_B`, `Sigma_A_last`, and `Sigma_B_last`
matrix entries for every method and lag. The displayed values are rounded to
four decimal places; the estimator calls retain their original MCMC settings.
