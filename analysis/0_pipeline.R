# PANDAnet0 - Script to run the pipeline for PANDAnet - Giulia, March 2024
# (Set working directory if necessary)
# ==== 1: Data setup ====
source("analysis/1_data_setup.R")

# ==== 2: Linear Mixed Models ====
source("analysis/2_lmm.R")

# ==== 3: Contemporaneous Networks ====
source("analysis/3_contemporaneous_networks.R")

# ==== 4: Temporally lagged networks ====
source("analysis/4_temporally_lagged_networks.R")

# ==== 5: Plots ====
source("analysis/5_plots.R")
