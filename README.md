The effect of sertraline on networks of mood and anxiety symptoms:
secondary analyses of the PANDA randomised controlled trial
================
Giulia Piazza

# Description

This repository contains code used in the analysis for the manuscript
“The effect of sertraline on networks of mood and anxiety symptoms:
secondary analyses of the PANDA randomised controlled trial”. The code
in the folder “analysis” can be run in R or R Studio. Scripts are
numbered, and should be run in order (this can be done by running
[pipeline.R](https://github.com/giuliapiazza18/PANDAnet-2/blob/main/analysis/0_pipeline.R)).
Tables, network matrices and figures are saved in the folder “results”.

## 0. Pipeline

This
[script](https://github.com/giuliapiazza18/PANDAnet-2/blob/main/analysis/0_pipeline.R)
runs all subsequent scripts in the analysis.

## 1. Data Setup

This
[script](https://github.com/giuliapiazza18/PANDAnet-2/blob/main/analysis/1_data_setup.R)
loads and inspects data, tidying data frames to be used in linear mixed
models and network analyses.

## 2. Linear Mixed Models

This
[script](https://github.com/giuliapiazza18/PANDAnet-2/blob/main/analysis/2_lmm.R)
estimates linear mixed models.

## 3. Contemporaneous Networks

This
[script](https://github.com/giuliapiazza18/PANDAnet/blob/main/analysis/3_contemporaneous_networks.R)
estimates contemporaneous networks and compares network structures with
the Network Comparison Test.

## 4. Temporally lagged Networks

This
[script](https://github.com/giuliapiazza18/PANDAnet-2/blob/main/analysis/4_temporally_lagged_networks.R)
estimates temporally lagged networks, and compares network structures
with model comparisons (fit indices, chi-square comparison).

## 5. Plots

This
[script](https://github.com/giuliapiazza18/PANDAnet-2/blob/main/analysis/5_plots.R)
plots and saves relevant figures.
