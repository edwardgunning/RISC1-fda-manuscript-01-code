Code for the paper ‘Analysing kinematic data from recreational runners
using functional data analysis’
================

## Repository Structure:

- :open_file_folder: **code**
  - :open_file_folder: **analysis** – scripts used to perform the data
    analysis.
    - :page_facing_up: [01 - Exploratory Plot for
      Introduction](code/analysis/BFMM-introduction-plot.R)
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
    - :page_facing_up: Scripts for comparison with existing approaches
      ([1](code/analysis/BFMM-multifamm-comparison.R),
      [2](code/analysis/BFMM-multifamm-comparison-figures.R),
      [3](code/analysis/BFMM-multifamm-comparison-boot-02.R),
      [4](code/analysis/BFMM-multifamm-comparison-boot-figures-02.R),
      [5](code/analysis/BFMM-fui-comparison.R),
      [6](code/analysis/BFMM-fui-comparison-figure.R))
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
    - :open_file_folder: [Functions for the FUI comparison (with a
      note)](code/functions/FUI-functions/)
    - :page_facing_up: [Function for multiFAMM
      comparison](code/functions/rough_fit_mfamm_model.R)
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

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Monterey 12.2.1
    ## 
    ## Matrix products: default
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_IE.UTF-8/en_IE.UTF-8/en_IE.UTF-8/C/en_IE.UTF-8/en_IE.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  splines   stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] lme4_1.1-30       fda_5.5.1         deSolve_1.30      fds_1.8          
    ##  [5] RCurl_1.98-1.6    rainbow_3.6       pcaPP_1.9-74      MASS_7.3-55      
    ##  [9] Matrix_1.4-0      data.table_1.14.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.10        pracma_2.3.8       pillar_1.7.0       compiler_4.1.2    
    ##  [5] nloptr_2.0.0       bitops_1.0-7       ks_1.13.4          tools_4.1.2       
    ##  [9] boot_1.3-28        mclust_5.4.9       tibble_3.1.6       lifecycle_1.0.3   
    ## [13] nlme_3.1-155       lattice_0.20-45    pkgconfig_2.0.3    rlang_1.0.6       
    ## [17] DBI_1.1.2          cli_3.6.0          rstudioapi_0.13    mvtnorm_1.1-3     
    ## [21] stringr_1.4.0      dplyr_1.0.8        cluster_2.1.2      generics_0.1.2    
    ## [25] vctrs_0.5.1        tidyselect_1.1.1   rprojroot_2.0.2    grid_4.1.2        
    ## [29] glue_1.6.2         hdrcde_3.4         here_1.0.1         R6_2.5.1          
    ## [33] fansi_1.0.2        minqa_1.2.4        purrr_0.3.4        magrittr_2.0.2    
    ## [37] ellipsis_0.3.2     matrixcalc_1.0-5   assertthat_0.2.1   colorspace_2.0-3  
    ## [41] utf8_1.2.2         KernSmooth_2.23-20 stringi_1.7.6      crayon_1.5.0

We have stored `.Random.seed` on each iteration of the simulation in a
list called `simulation_seeds`. The seed can be set to produce any
iteration by using the command (described
[here](https://stackoverflow.com/questions/19614314/can-i-get-seed-somehow)):

``` r
# to set seed as it was on iteration i
given_seed <- simulation_seeds[[i]]
.Random.seed <- given_seed
```
