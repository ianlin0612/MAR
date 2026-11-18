# Matrix autoregression estimators

This repository provides estimators for matrix autoregressive (MAR) models and
code for three simulation experiments in the paper.

## Core functions

The complete estimator code is in [`R/mar_estimators.R`](R/mar_estimators.R).
It defines two functions that take a time series of matrices and estimate the
row and column coefficient matrices of a MAR model. The functions also return
posterior inclusion probabilities (PIPs), which indicate how often each
coefficient was included in the sampled models, and estimates of the row and
column error covariance matrices.

Run `source("R/mar_estimators.R")` in R to load these functions before calling
them. For package installation, the meaning of each input, how to call the
functions, and how to read their outputs, see
[`docs/MAR_ESTIMATORS.md`](docs/MAR_ESTIMATORS.md). This document explains the
functions; `R/mar_estimators.R` is the file containing their executable code.

## Simulation examples

Each example corresponds to one experiment in the paper. It contains the code
for generating the true parameters and simulated data, the inputs passed to
the core functions, and the results from one simulation run for comparison.

| Paper experiment | Example file |
| --- | --- |
| Simulation 1 | [`docs/EXAMPLE_SIMULATION1.md`](docs/EXAMPLE_SIMULATION1.md) |
| Simulation 2 | [`docs/EXAMPLE_SIMULATION2.md`](docs/EXAMPLE_SIMULATION2.md) |
| Supplement S2 | [`docs/EXAMPLE_SUPPLEMENT_S2.md`](docs/EXAMPLE_SUPPLEMENT_S2.md) |

## Running an example

Install the packages listed in [`docs/MAR_ESTIMATORS.md`](docs/MAR_ESTIMATORS.md#packages),
then start a fresh R session.
Run the following in order:

1. Load the core-function code by running `source("R/mar_estimators.R")`.
2. Open the example file from the table above or the
   [Examples section](docs/MAR_ESTIMATORS.md#examples) of `MAR_ESTIMATORS.md`,
   then run its **Data generation** code block. This sets the true parameters
   and generates the training data.
3. Run the **Core-function inputs** code block in the same example. It passes
   the generated data and experiment inputs to the core functions and stores
   the fits in `fits_by_rep`.

Keep the same R session for steps 1 through 3 so that the functions and generated data
remain available when the estimator calls run.

Each example also provides a **One-seed core-function output** section with
coefficient estimates, PIPs, and covariance estimates from one simulation run
for reference.
