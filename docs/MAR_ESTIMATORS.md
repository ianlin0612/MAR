# Core MAR estimators

The complete core-function code is in
[`R/mar_estimators.R`](../R/mar_estimators.R).

## Contents

- [Description](#description): Model and estimators
- [Packages](#packages): Required R packages
- [Usage](#usage): Function call signatures
- [Inputs](#inputs): Function arguments
- [Output](#output): Returned values
- [Examples](#examples): Simulation experiment files

## Description

The functions below estimate the row coefficient matrices $A_{\ell}$ and
column coefficient matrices $B_{\ell}$ at each lag in the matrix autoregression

```math
Y_t = \sum_{\ell=1}^{P} A_{\ell} Y_{t-\ell} B_{\ell}^{\mathsf{T}} + E_t.
```

`fit_mar_pointmass()` uses a point-mass spike-and-slab prior. With individual
coefficient groups it reports the CGS method; with grouped coefficients it
reports GGS. `fit_mar_ssvs()` uses a continuous spike-and-slab prior and
reports SSVS. Both functions estimate separable row and column innovation
covariances and return posterior inclusion probabilities (PIPs) and sparse
coefficient estimates for every lag.

## Packages

Required R packages: `MASS`, `MCMCpack`, and `Matrix`.

If needed, install them first:

```r
install.packages(c("MASS", "MCMCpack", "Matrix"))
```

```r
library(MASS)
library(MCMCpack)
library(Matrix)
```

## Usage

### `fit_mar_pointmass`

```r
fit_mar_pointmass(Y, P = 1L, partition_A = NULL, partition_B = NULL,
    A_init = NULL, B_init = NULL, Sigma_A_init = NULL, Sigma_B_init = NULL,
    n_iter = 10000L, burnin = 5000L, warm_iter = 20L,
    tau_A = 4, tau_B = 4, prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2), nu_A = NULL, nu_B = NULL,
    scale_init = NULL, scale_df = NULL, eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = c(">", ">="), rescale = TRUE,
    random_lag_scan = FALSE, init_candidates = NULL,
    singleton_update = c("auto", "scalar", "group_1x1"),
    seed = NULL, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL,
    mcmc_clock = NULL)
```

### `fit_mar_ssvs`

```r
fit_mar_ssvs(Y, P = 1L,
    A_init = NULL, B_init = NULL, Sigma_A_init = NULL, Sigma_B_init = NULL,
    n_iter = 2250L, burnin = 750L, warm_iter = 10L,
    init_method = c("map", "gibbs", "none"),
    tau0_A = 0.05, tau1_A = 2, tau0_B = tau0_A, tau1_B = tau1_A,
    prior_A = list(diag = 0.4, off = 0.2), prior_B = list(diag = 0.4, off = 0.2),
    fixed_theta = TRUE, theta_A_init = NULL, theta_B_init = NULL,
    alpha_theta = 1, beta_theta = 1, warm_fixed_theta = TRUE,
    nu_A = NULL, nu_B = NULL, scale_init = 1, scale_df = NULL, eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = c(">", ">="), rescale = TRUE,
    scale_update = c("before_rescale", "after_rescale"),
    draw_initial_indicators = TRUE, seed = NULL, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL,
    mcmc_clock = NULL)
```

## Inputs

### Inputs shared by both functions

| Input | Definition |
| --- | --- |
| `Y` | Finite numeric training array of dimensions `G x N x T_train`; `T_train > P`. Each time slice is one `G x N` observation. |
| `P` | Positive integer autoregressive order. |
| `A_init`, `B_init` | Initial coefficient matrices of size `G x G` and `N x N`. Supply one matrix for every lag, or a list of `P` matrices. `NULL` uses identity matrices. |
| `Sigma_A_init`, `Sigma_B_init` | Initial symmetric positive-definite row and column covariance matrices, of sizes `G x G` and `N x N`; `NULL` uses identities. |
| `n_iter` | Total number of main MCMC sweeps, including burn-in. |
| `burnin` | Initial main-chain sweeps discarded from posterior summaries; must be smaller than `n_iter`. |
| `warm_iter` | Number of initialization sweeps before the main MCMC chain. For point-mass estimation, this applies to each initial-state candidate; for SSVS, the warm-up method is selected by `init_method`. |
| `prior_A`, `prior_B` | Prior inclusion probabilities for row and column coefficients. Each is a scalar or `list(diag = ..., off = ...)` with values strictly between 0 and 1. |
| `nu_A`, `nu_B` | Inverse-Wishart degrees of freedom for row and column covariances. `NULL` means `G + 2` and `N + 2`; each must exceed its dimension minus one. |
| `scale_init` | Positive initial shared covariance scale. The default is `(G + N) / 2` for point-mass estimation and `1` for SSVS. |
| `scale_df` | Positive numeric vector of length two used by the shared-scale update; `NULL` means `c(nu_A, nu_B)`. |
| `eta_shape`, `eta_rate` | Positive shape and rate parameters of the gamma prior on the shared scale. |
| `pip_cut` | Threshold in `[0, 1]` for selecting coefficients from PIPs in `A_hat` and `B_hat`. |
| `pip_rule` | `">"` or `">="`; defines whether a PIP must exceed or may equal `pip_cut`. |
| `rescale` | Logical flag: balance norms of paired coefficient and covariance factors after each sweep. |
| `seed` | Optional nonnegative integer RNG seed before initialization. |
| `warm_seed` | Optional nonnegative integer RNG seed for warm-up. |
| `sampling_seed` | Optional nonnegative integer RNG seed reset before the main MCMC chain. |
| `keep_draws` | Logical flag: include retained MCMC draws in the returned list. |
| `verbose_every` | Positive integer progress-message interval in sweeps, or `NULL` for no messages. |
| `mcmc_clock` | Optional clock object with `start()`, `stop()`, and `cancel()` methods, timed around the main MCMC loop. |

### Additional point-mass inputs

| Input | Definition |
| --- | --- |
| `partition_A`, `partition_B` | Lists of coefficient index groups partitioning `1:G` and `1:N`, with every index used once. `NULL` uses singleton groups. Supplying groups changes the reported method from CGS to GGS. |
| `tau_A`, `tau_B` | Positive slab standard deviations for row and column coefficients. Each is a scalar or `list(diag = ..., off = ...)`. |
| `random_lag_scan` | Logical flag: randomly order lag updates in each sweep when `TRUE`. |
| `init_candidates` | Optional nonempty list of candidate initial states, each containing `A` and `B` matrices (one or `P` per side). Each receives `warm_iter` sweeps; the candidate with the smallest residual sum of squares is selected. `NULL` uses only `A_init` and `B_init`. |
| `singleton_update` | `"scalar"`, `"group_1x1"`, or `"auto"`. Auto uses scalar updates for a `NULL` partition and group updates otherwise. |


### Additional SSVS inputs

| Input | Definition |
| --- | --- |
| `init_method` | `"map"`, `"gibbs"`, or `"none"` warm-up method. `"none"` skips warm-up even when `warm_iter` is positive. |
| `tau0_A`, `tau0_B` | Positive spike standard deviations for row and column coefficients. |
| `tau1_A`, `tau1_B` | Positive slab standard deviations, each greater than its corresponding spike scale. |
| `fixed_theta` | Logical flag: keep inclusion probabilities fixed during MCMC. |
| `theta_A_init`, `theta_B_init` | Optional initial inclusion probabilities: scalar, `d x d` matrix, or `d x (d * P)` matrix for `d = G` or `N`. `NULL` uses `prior_A`/`prior_B` when `fixed_theta = TRUE`, otherwise the beta-prior mean. |
| `alpha_theta`, `beta_theta` | Positive beta-prior shape parameters for inclusion probabilities when `fixed_theta = FALSE`. |
| `warm_fixed_theta` | Logical flag: hold inclusion probabilities fixed during warm-up when `fixed_theta = FALSE`. |
| `scale_update` | `"before_rescale"` or `"after_rescale"`; selects when the shared scale is updated relative to factor balancing. |
| `draw_initial_indicators` | Logical flag: sample initial inclusion indicators before warm-up. |


## Output

Each function returns a list.

| Returned field | Definition |
| --- | --- |
| `method`, `P` | Method label (`CGS`, `GGS`, or `SSVS`) and lag order. |
| `A_hat`, `B_hat` | Named lists of `P` PIP-thresholded coefficient estimates; each row matrix is `G x G`, each column matrix `N x N`. |
| `PIP_A`, `PIP_B` | Named lists of matching posterior inclusion probability matrices. Point-mass PIPs use nonzero draws; SSVS PIPs use sampled inclusion indicators. |
| `Sigma_A_last`, `Sigma_B_last` | Covariance-factor matrices from the final sweep. |
| `scale_last` | Shared covariance scale from the final sweep. |
| `runtime_sec` | Elapsed seconds for the main MCMC loop. |
| `settings` | List of selected settings recorded by the estimator. It is not a complete copy of all call inputs. |
| `A_chain`, `B_chain` | Only if `keep_draws = TRUE`: numeric arrays of dimensions `(n_iter - burnin) x P x G x G` and `(n_iter - burnin) x P x N x N`. |
| `scale_chain` | Only if `keep_draws = TRUE`: retained shared-scale draws. |

## Examples

The three example documents contain the complete data-generation and
estimator-call blocks. In a fresh R session at the repository root, run
`source("R/mar_estimators.R")` first. Then run each example document's
data-generation block followed by its estimator-call block.
The second code block stores all estimator outputs in `fits_by_rep`; the
one-seed output section shows how to print the first replication. Set `n_reps`
at the top of the first code block to choose the number of replications.

| Experiment | Example document | Methods |
| --- | --- | --- |
| Simulation 1, MAR(1), `20 x 20 x 90` | [`EXAMPLE_SIMULATION1.md`](EXAMPLE_SIMULATION1.md) | CGS, SSVS |
| Simulation 2, MAR(3), `15 x 10 x 450` | [`EXAMPLE_SIMULATION2.md`](EXAMPLE_SIMULATION2.md) | GGS, SSVS |
| Supplement S2, MAR(1), `24 x 28 x 160` | [`EXAMPLE_SUPPLEMENT_S2.md`](EXAMPLE_SUPPLEMENT_S2.md) | CGS, GGS |
