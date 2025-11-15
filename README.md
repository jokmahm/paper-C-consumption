# Functional Response Model for Empirical Cod Consumption

This repository contains the code used in **Paper C of my PhD thesis**, titled:

**"A Functional Response Model Fit to Empirical Consumption Data"**

The work develops and fits an analytical functional response model to empirical cod consumption data using **R** and **TMB (Template Model Builder)**. The estimation uses C++ templates for the likelihood, R scripts for running and analysing the model, and optional MATLAB scripts for generating visualisations.

All code required to reproduce the results in the paper is included in this repository.

---

## How to Run the R + TMB Code

1. Install required R packages:
    ```r
    install.packages(c("TMB", "ggplot2", "dplyr"))
    ```

2. Compile the TMB model:
    ```r
    library(TMB)
    compile("cpp/Consumption.cpp")
    dyn.load(dynlib("cpp/Consumption"))
    ```

3. Run the main estimation script:
    ```r
    source("R/1_run_model_no_cod_montecarlo.R")
    ```

4. Run alternative model versions:
    - Selective temperature:
        ```r
        source("R/1_run_model_no_cod_montecarlo_2_select.R")
        ```
    - Selective model with epsilon:
        ```r
        source("R/1_run_model_no_cod_montecarlo_2_select_eps.R")
        ```

---

## MATLAB Figures (Optional)

To generate the MATLAB visualisations used in the paper:
```matlab
cd matlab/3D-surface
run("main.m")
