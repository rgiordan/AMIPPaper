#!/usr/bin/env Rscript
################################################################################
# Process the output of the scripts
# microcredit_vb.py and sensitivity.py into a more convenient Rdata format.

# microcredit tailored hierarchical pdf for profit with lognormal tails
# Rachael Meager
# April 2016 

library(optparse)

if (interactive()) {
  git_repo_loc <- system("git rev-parse --show-toplevel", intern=TRUE)
  num_advi_draws <- 20
  python_path <- "python"
} else {
  option_list <- list(
    make_option("--git_repo_loc", type="character", default=NULL, 
                help="Paper github repo", metavar="character"),
    make_option("--python_path", type="character", default=NULL, 
                help="Path to virtual environment python", metavar="character"),
    make_option("--num_draws", type="integer", default=NULL, 
                help="Number of draws", metavar="integer")
  )
  opt_parser <- OptionParser(option_list=option_list);
  opt <- parse_args(opt_parser);
  git_repo_loc <- opt$git_repo_loc
  num_advi_draws <- opt$num_draws
  python_path <- opt$python_path
}
if (!dir.exists(git_repo_loc)) {
  stop(sprintf("Repository directory %s does not exist.", git_repo_loc))
}

# Set this to your path.
base_dir <- file.path(git_repo_loc, "examples/microcredit_mixture")

# Set these to the values for your particular analysis.
use_smuber <- "False"

print(sprintf("Using python %s", python_path))

# These two sanity checks are not essential for the final paper
# and are time-consuming.  For now let's remove them for the
# purpose of reproducibility.
save_stan_comparison <- FALSE
save_freq_check <- FALSE


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

setwd(file.path(base_dir, "python"))
source("advi_sensitivity_lib.R")

# Load the raw data
data_filename <- file.path(base_dir, "R/data/microcredit_project_data_cleaned.csv")
data <- read.csv(data_filename, header=TRUE)


initial_fit_filename <- file.path(
    base_dir, "python/output/",
    sprintf("microcredit_project_advi_%ddraws_smuber%s.npz", 
            num_advi_draws, use_smuber))

sens_filename <- sprintf(
    "output/microcredit_project_weight_sensitivity_%ddraws_smuber%s.npz",
    num_advi_draws, use_smuber)

#python_path <- file.path(base_dir, "venv/bin/python3")
py_main <- InitializePythonAndLoadFit(
  initial_fit_filename, sens_filename, base_dir, python_path=python_path)
DefineCustomParameterFunctions(py_main)

# Load the Stan results
if (save_stan_comparison) {
  stan_output <- new.env()
  load(file.path(base_dir, "R/output/",
                 "microcredit_profit_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData"),
       env=stan_output)
  sampling_result <- stan_output$stan_fit
  stan_summary <- summary(sampling_result)$summary
}


############################
# Make an R list for saving results.

save_list <- list()
save_filename <- file.path(
    base_dir, "data/", 
    sprintf("microcredit_vb_%s_combined_results.Rdata", analysis_desc))

##########################################
# Save the raw data.  (Nice to have it compressed for faster loading rather than in a csv)

save_list[["data"]] <- data
save_list[["initial_fit_filename"]] <- initial_fit_filename
save_list[["sens_filename"]] <- sens_filename
save_list[["data_filename"]] <- data_filename

################################
# Compare MCMC and ADVI means

GatherADVIResult <- function(result) {
    GatherMicrocreditResult(result,
                            param_names=py_main$flat_names,
                            increment_indices=TRUE) %>%
        mutate(method="advi")
}

advi_df <- bind_rows(
    GatherADVIResult(py_main$mean_params) %>% mutate(metric="mean"),
    GatherADVIResult(sqrt(diag(py_main$lr_cov))) %>% mutate(metric="sd"))

if (save_stan_comparison) {
  stan_df <- bind_rows(
    GatherMicrocreditResult(t(stan_summary[, "mean"])) %>%
      mutate(metric="mean", method="mcmc"),
    GatherMicrocreditResult(t(stan_summary[, "sd"])) %>%
      mutate(metric="sd", method="mcmc"))
  
  check_df <-
    inner_join(advi_df, stan_df, by=c("param", "k", "p", "m", "s", "metric"))
  stopifnot(nrow(check_df) == nrow(advi_df))
  
  
  comp_df <-
    bind_rows(advi_df, stan_df) %>%
    mutate(method_metric=paste(method, metric, sep="_")) %>%
    dcast(param + k + p + m + s ~ method_metric, value.var="value")
  head(comp_df)
  
  save_list[["comp_df"]] <- comp_df
  save_list[["num_mcmc_draws"]] <- ncol(py_main$base_draws)
}


# Save some timing results

# The total amount of time in seconds to optimize the ADVI objective.
# As of writing I had not save the amount of time to compute the Hessian, but
# it was relatively small, about a minute.
save_list[["advi_time"]] <- py_main$advi_time


########################################################
# Estimate error due to using only a few ADVI samples

if (save_freq_check) {
  # Load the frequentist sensitivity results
  py_main$freq_filename <- sprintf(
    "output/microcredit_project_advi_frequentist_num_draw_check_%ddraws_smuber%s.npz",
    num_advi_draws, use_smuber)
  
  print(py_main$sens_filename)
  reticulate::py_run_string("
freq_dict = np.load(freq_filename)
rel_error_flat = freq_dict['rel_error_flat']
rel_freq_error = param_pattern.fold(rel_error_flat, free=False)
")
  
  save_list[["relative_freq_error"]]  <-
    GatherADVIResult(py_main$rel_error_flat) %>%
    mutate(metric="advi_se_over_post_sd")
}



#################################
# Influence

dim(py_main$mean_influence_mat)
mean_influence_mat <- as.matrix(t(py_main$mean_influence_mat))
colnames(mean_influence_mat) <- py_main$flat_names

save_list[["mean_influence_mat"]] <- mean_influence_mat
save_list[["mean_params"]] <- py_main$mean_params
save_list[["lr_cov"]] <- py_main$lr_cov
save_list[["param_names"]] <- py_main$flat_names


###################
# Get the derivative of custom functions, in particular of log(sd_tau).

custom_results <- EvalCustomParameters(py_main, py_main$advi_params)

# Sanity check --- these should match approximately
sd_tau_inds <- py_main$flat_names %in% sprintf("sd_tau[%d]", c(0, 1))
log(py_main$mean_params[sd_tau_inds])
custom_results$mean

save_list[["custom_param_names"]] <- c("log_sd_tau[0]", "log_sd_tau[1]")
save_list[["custom_influence_mat"]] <- t(custom_results$infl)
save_list[["custom_advi_mean"]] <- custom_results$mean
save_list[["custom_advi_cov"]] <- custom_results$sd


###############################
# Save extra information for refitting.

save_list[["advi_params_free"]] <- py_main$advi_params_free
save_list[["base_draws"]] <- py_main$base_draws


#############################
# Save.

print(sprintf("Saving to %s", save_filename))
save(save_list, file=save_filename)

print("Success!")