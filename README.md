Code for the paper ‘Analysing kinematic data from recreational runners
using functional data analysis’
================

## Repository Structure:

- :open_file_folder: **code**
  - :open_file_folder: **analysis** – scripts used to perform the data
    analysis.
    - :page_facing_up: [01 - Exploratory Plot for
      Introduction](code/analysis/BFMM-paper-basis-transformation.R)
    - :page_facing_up: [02 - Basis
      Transformation](code/analysis/BFMM-paper-basis-transformation.R)
    - :page_facing_up: [03 - Model
      Fitting](code/analysis/BFMM-paper-modelling.R)
    - :page_facing_up: [04 - Bootstrap and Wald
      Inference](code/analysis/BFMM-paper-bootstrap.R)
    - :page_facing_up: [05 - Plotting Fixed-Effects
      Results](code/analysis/BFMM-paper-fixef-results.R)
    - :page_facing_up: [06 - Plotting Random-Effects
      Results](code/analysis/BFMM-paper-covariance-results.R)
    - :page_facing_up: [07 - Additional Exploratory Analysis of
      Random-Effects
      Results](code/analysis/BFMM-paper-covariance-extra.R)
    - :page_facing_up: [08 - Analysis of
      ICC](code/analysis/BFMM-icc-analysis.R)
  - :open_file_folder: **functions** – custom functions used to perform
    the data analysis and simulation
    - :page_facing_up: [Custom `ggplot2` theme for
      figures](code/functions/theme_gunning.R)
    - :page_facing_up: [Functions to help storage and manipulation of
      `fda::fd` objects](code/functions/functions-helper-smoothing.R)
    - :page_facing_up: [Function to calculate uncentered FPCA
      scores](code/functions/function-project-mean-onto-fpcs.R)
    - :page_facing_up: [Function to calculate unstructured covariance
      estimates (based off the `denseFLMM`
      package)](code/functions/functions-unstructured-covariance.R)
    - :page_facing_up: :page_facing_up: :page_facing_up:
      :page_facing_up: Tests for the unstructured covariance estimator
      [(1)](code/functions/cov_unstruct_test-01.R)
      [(2)](code/functions/cov_unstruct_test-02.R)
      [(3)](code/functions/cov_unstruct_test-03.R)
      [(4)](code/functions/cov_unstruct_test-04.R)
    - :page_facing_up: [Function to extract an estimated heteroscedastic
      residual covariance matrix from an `nlme`
      object](code/functions/function-get-residual-covariance-matrix.R)
    - :page_facing_up: [Function to calculate Monte Carlo Standard
      Errors (SEs) for Coverage](code/functions/binomial_se.R)
    - :page_facing_up: [Tests Monte Carlo SEs for
      Coverage](code/functions/binomial_se_tests.R)
  - :open_file_folder: **simulation** – scripts to generate and perform
    the short simulation in the paper
    - :page_facing_up: [01 - Fit a simple model and generate empirical
      parameters](code/simulation/BFMM-paper-get-simulation-parameters.R)
    - :page_facing_up: [02 - Generate Eigenfunctions for the
      Simulation](code/simulation/BFMM-paper-generate-efuns-simulation.R)
    - :page_facing_up: [03 - Functions to Generate
      Data](code/simulation/BFMM-paper-generate-simulated-data.R)
    - :page_facing_up: [04 - Code to Run and Save
      Simulation](code/simulation/BFMM-paper-tidied-simulation.R)
    - :page_facing_up: [05 - Assess Estimation and Plot
      Results](code/simulation/BFMM-paper-simulation-result-plot.R)
    - :page_facing_up: [06 - Assess Coverage and Tabulate
      Results](code/simulation/BFMM-paper-simulation-coverage-tables.R)

## Timing Results:

It took **15.32 minutes** to run the bootstrap analysis using **2500
bootstrap replicates** with computing shared across **8 cores** on a
2019 MacBook Pro with a 2.4 GHz Quad-Core Intel Core i5 processer and 8
GB of memory (code ran on 2023-02-17 13:42:16).

## Reproducibility

We have stored `.Random.seed` on each iteration of the simulation in a
list called `simulation_seeds`. The seed can be set to produce any
iteration by using the command (described
[here](https://stackoverflow.com/questions/19614314/can-i-get-seed-somehow)):

``` r
# to set seed as it was on iteration i
given_seed <- simulation_seeds[[i]]
.Random.seed <- given_seed
```
