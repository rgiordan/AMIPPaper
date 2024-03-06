
#######################################################
# Compare the output of VB to MCMC using the output of process_result_for_paper.R
# I can't find where I did this before, so this might be redundant...

# Set this to your path.
base_dir <- "/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench/examples/microcredit_mixture/"

# Set these to the values for your particular analysis.
use_smuber <- "False"
num_advi_draws <- 20


##########################
library(ggplot2)
library(rstan)
library(reticulate)
library(dplyr)
rstan_options(auto_write = TRUE)
cores <- as.numeric(1)
library(reshape2)
library(tidybayes)

library(tidyr)
library(zaminfluence)
library(stringr)

analysis_desc <- sprintf("draws%d_smuber%s", num_advi_draws, use_smuber)

save_filename <- file.path(base_dir, "data/", sprintf("microcredit_vb_%s_combined_results.Rdata", analysis_desc))
load(save_filename)

comp_df <- save_list$comp_df

head(comp_df)

# Compare the means
ggplot(comp_df) +
    geom_point(aes(y=advi_mean, x=mcmc_mean))

# Comapre the SDs.  Note that "ADVI sd" actually means LRVB sd.
ggplot(comp_df) +
    geom_point(aes(y=advi_sd, x=mcmc_sd)) +
    geom_abline(aes(slope=1, intercept=0))


# Compare the sampling error realtive to the posterior sd.
relative_freq_error <- save_list$relative_freq_error
ggplot(relative_freq_error) +
    geom_histogram(aes(x=value), bins=50)