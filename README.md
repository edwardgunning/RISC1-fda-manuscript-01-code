Code for the paper ‘Analysing kinematic data from recreational runners
using functional data analysis’
================

## Repository Structure:

- :open_file_folder: **code**
  - :open_file_folder: **analysis** – scripts used to perform the data
    analysis.
    - :page_facing_up: [01 - Basis
      Transformation](code/analysis/BFMM-paper-basis-transformation.R)
    - :page_facing_up: [02 - Model
      Fitting](code/analysis/BFMM-paper-modelling.R)
    - :page_facing_up: [03 - Bootstrap and Wald
      Inference](code/analysis/BFMM-paper-bootstrap.R)
    - :page_facing_up: [04 - Plotting Fixed-Effects
      Results](code/analysis/BFMM-paper-bootstrap.R)
    - :page_facing_up: [05 - Plotting Random-Effects
      Results](code/analysis/BFMM-paper-covariance-results.R)
  - :open_file_folder: **functions** – custom functions used to perform
    the data analysis.
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

## Timing Results:

It took **15.7 minutes** to run the bootstrap analysis using **2500
bootstrap replicates** with computing shared across **8 cores** on a
2019 MacBook Pro with a 2.4 GHz Quad-Core Intel Core i5 processer and 8
GB of memory (code ran on 2023-02-07 15:54:28).
