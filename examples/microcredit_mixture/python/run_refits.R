#!/usr/bin/env Rscript

################################################################################
# Define, print commands for, and summarize refits using output of the notebooks
# microcredit_vb.ipynb and sensitivity.ipynb, and the R script save_initial_fit_as_rdata.R.


library(optparse)

if (interactive()) {
  git_repo_loc <- system("git rev-parse --show-toplevel", intern=TRUE)
  num_advi_draws <- 30
  python_path <- "/home/rgiordan/miniforge3/envs/amip-paper-2023/bin/python"
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


# Set this to your path.
base_dir <- file.path(git_repo_loc, "examples/microcredit_mixture/")

# Libraries
library(ggplot2)
library(tidyverse)
library(zaminfluence)

setwd(file.path(base_dir, "python"))
source(file.path(base_dir, "python/advi_sensitivity_lib.R"))


########################
# Load the data.

# Set these to the values for your particular analysis.

# Whether or not to use the "smuber" loss.  I think we're going with False for the paper.
use_smuber <- "False"


# The filename with the data we're going to use.  This was produced by process_results_for_paper.
analysis_desc <- sprintf("draws%d_smuber%s", num_advi_draws, use_smuber)
save_filename <- file.path(
    base_dir, "data/",
    sprintf("microcredit_vb_%s_combined_results.Rdata", analysis_desc))

# Load the analysis.
load(save_filename)

# Load the Python results.  Unfortunately, as it's currently set up,
# we need the Python results to compute the custom quantities,
# which seems unnecessarily complicated and not ideal.
py_main <- InitializePythonAndLoadFit(
    save_list[["initial_fit_filename"]], 
    save_list[["sens_filename"]],
    base_dir, 
    python_path=python_path)
DefineCustomParameterFunctions(py_main)

# Define an output file
output_dir <- file.path(git_repo_loc, "writing/output/applications_data/microcredit_mixture")
if (!dir.exists(output_dir)) {
  stop(sprintf("Output directory %s does not exist.", output_dir))
}

##########################################
# Convert to ModelGrads object.

param_names=save_list[["param_names"]]
mean_params=save_list[["mean_params"]]
mean_influence_mat=save_list[["mean_influence_mat"]]
lr_cov=save_list[["lr_cov"]]
custom_param_names=save_list[["custom_param_names"]]
custom_advi_mean=save_list[["custom_advi_mean"]]
custom_influence_mat=save_list[["custom_influence_mat"]]
custom_advi_cov=save_list[["custom_advi_cov"]]









model_grads <- PythonToModelGrads(
    param_names=save_list[["param_names"]],
    mean_params=save_list[["mean_params"]],
    mean_influence_mat=save_list[["mean_influence_mat"]],
    lr_cov=save_list[["lr_cov"]],
    custom_param_names=save_list[["custom_param_names"]],
    custom_advi_mean=save_list[["custom_advi_mean"]],
    custom_influence_mat=save_list[["custom_influence_mat"]],
    custom_advi_cov=save_list[["custom_advi_cov"]])


# Define some interesting parameters.  Note that the parameters are zero-indexed here.
param_names <- c(
    sprintf("tau[%d]", 0:1),
    sprintf("log_sd_tau[%d]", 0:1),
    as.character(outer(0:6, 0:1, function(i,j) sprintf("tau_k[%d,%d]", i,j))))

nrow(model_grads$beta_grad)
length(model_grads$weights)

options(warn=2)
class(model_grads)

for (param in param_names) {
    model_grads <-
        model_grads %>%
        AppendTargetRegressorInfluence(param)
}


# TODO: update this one too
signals <- GetInferenceSignals(model_grads)


# signals_list <- list()
# for (param in sprintf("tau[%d]", 0:1)) {
#     signals_list[[param]] <-
#         GetInferenceSignals(model_grads$param_infl_list[[param]])
# }



############################################################################
# Save some weights for re-fitting.


SaveWeights <- function(weights, weight_filename) {
    py_main$weight_filename <- weight_filename
    print(sprintf("Saving to %s", py_main$weight_filename))
    py_main$w <- weights
    
    reticulate::py_run_string("
save_dict = {}
save_dict['weights'] = w
save_dict['out_filename'] = weight_filename
save_dict['advi_params_free'] = advi_params_free
save_dict['base_draws'] = base_draws
np.savez_compressed(**save_dict, file=weight_filename)
")
}


GetRefitFilename <- function(param, change, num_draws=num_advi_draws) {
    output_filename <- file.path(
        base_dir, "python/reweightings",
        sprintf("refit_%s_%s_%ddraws.npz", param, change, num_draws))
    return(output_filename)    
}


GetRefitCommand <- function(weight_filename, output_filename) {
    script_path <- file.path(base_dir, "python")
    cmd <- paste0(
        python_path, " \\\n",
        script_path, "/refit.py \\\n",
        "--initial_fit_filename=", save_list[["initial_fit_filename"]], " \\\n",
        "--data_filename=", save_list[["data_filename"]], " \\\n",
        "--reweight_filename=", weight_filename, " \\\n",
        "--output_filename=", output_filename)
    return(cmd)
}


#############################################################
# Generate weights and commands for running the re-fits


amis_list <- list()


# Refits for tau
for (param in sprintf("tau[%d]", 0:1)) {
    amis_list[[param]] <- list()
    for (change in c("sign", "sig", "both")) {
        apip <- signals[[param]][[change]]$apip
        if (!is.na(apip$n)) {
            amis <- apip$inds
            w_int <-
                GetWeightVector(
                    drop_inds=amis, num_obs=model_grads$model_fit$num_obs, bool=FALSE) %>%
                as.integer()
            stopifnot(sum(1 - w_int) == apip$n)
            weight_filename <-
                file.path(
                    base_dir, "python/reweightings",
                    sprintf("amis_weight_draws%d_%s_%s_smuber%s.npz",
                            num_advi_draws, param, change, use_smuber))
            cat(sprintf("Dropping %d for %s, %s", apip$n, param, change), "\n")
            SaveWeights(w_int, weight_filename)
            refit_filename <- GetRefitFilename(param, change)
            command <- GetRefitCommand(weight_filename, refit_filename)
            amis_list[[param]][[change]] <-
                list(amis=amis, n_drop=apip$n, prop_drop=apip$prop,
                     refit_filename=refit_filename,
                     command=command)
        } else {
            sprintf("APIP = NA for %s, %s", param, change)
        }
    }
}


names(model_grads$param_infls[[param]])

# Refits for log_tau_sd
# Drop up to 0.5% of points.
max_n_drop <- ceiling(0.005 * model_grads$model_fit$num_obs)
#n_drops <- ceiling(seq(0, max_n_drop, length.out=6))
n_drops <- c(0, 5, 10, 36, 71, 107, 142, 177)
stopifnot(max(n_drops) <= max_n_drop)
#prop_drops <- c(0.00, 0.001, 0.002, 0.003, 0.004, 0.005)
for (param in sprintf("log_sd_tau[%d]", 0:1)) {
    amis_list[[param]] <- list()
    for (n_drop in n_drops) {
        prop_drop <- n_drop / model_grads$model_fit$num_obs
        for (change_sign in c("neg", "pos")) {
            amis <- GetAMIS(
                qoi=model_grads$param_infls[[param]][["param"]], 
                sign=change_sign,
                n_drop=n_drop)
            if (!is.null(amis)) {
              w_int <-
                GetWeightVector(
                  drop_inds=amis, num_obs=model_grads$model_fit$num_obs, bool=FALSE) %>%
                as.integer()
              stopifnot(sum(1 - w_int) == n_drop)
              change <- sprintf("%s%d", change_sign, n_drop)
              weight_filename <-
                file.path(
                  base_dir, "python/reweightings",
                  sprintf("amis_weight_%s_%s_smuber%s.npz",
                          param, change, use_smuber))
              cat(sprintf("Dropping %d for %s, %s", apip$n, param, change), "\n")
              SaveWeights(w_int, weight_filename)
              refit_filename <- GetRefitFilename(param, change)
              command <- GetRefitCommand(weight_filename, refit_filename)
              amis_list[[param]][[change]] <-
                list(param=param, change=change,
                     amis=amis, n_drop=n_drop, prop_drop=prop_drop,
                     refit_filename=refit_filename,
                     command=command)
            }
        }
    }
}        

names(amis_list)
# Execute (or print) commands for missing refits.
execute <- TRUE
#for (param in sprintf("log_sd_tau[%d]", 0:1)) {
for (param in names(amis_list)) {
    amis_changes <- amis_list[[param]]
    for (change in names(amis_changes)) {
      print(sprintf("%s %s", param, change))
        amis <- amis_changes[[change]]
        if (!file.exists(amis$refit_filename)) {
            #cat("\n", param, change, "\n")
          cat("\n Executing\n", amis$command, " & \n\n")
          if (execute) {
            system(amis$command)
          } 
        }
    }
}


#############################################################
# Load the results of the refits

# A helper function to load the refits
LoadRefit <- function(refit_filename) {
    refit_filename <- as.character(refit_filename)
    if (!file.exists(refit_filename)) {
        cat("\nFile ", refit_filename, " does not exist; skipping.\n")
        return(NULL)
    }
    py_main$refit_filename <- refit_filename
    reticulate::py_run_string("
refit_dict = np.load(refit_filename)
refit_advi_params_free = refit_dict['advi_params_free']
refit_advi_params = advi_param_pattern.fold(refit_advi_params_free, free=True)
refit_mean_params = mm_lib.get_advi_mean(refit_advi_params, base_draws, param_pattern)
refit_lr_cov = refit_dict['lr_cov']
")
    custom_results <- EvalCustomParameters(py_main, py_main$refit_advi_params)
    # Need to add the custom parameters
    refit <- list(
        betahat=c(py_main$refit_mean_params, custom_results$mean),
        se=c(sqrt(diag(py_main$refit_lr_cov)), custom_results$sd))
    return(refit)
}


# Load refits from disk
if (!exists("refits_list")) {
    refits_list <- list()
}
for (amis_param in amis_list) {
    for (change_list in amis_param) {
        if (is.null(refits_list[[change_list$refit_filename]])) {
            cat("Loading ", change_list$refit_filename, "\n")
            refits_list[[change_list$refit_filename]] <-
                LoadRefit(change_list$refit_filename)
        }
    }
}


#############################################################
# Postprocess into a summary dataframe and save

refits_df <- data.frame()
for (param in names(amis_list)) {
    amis_param <- amis_list[[param]]
    for (change in names(amis_param)) {
        change_list <- amis_param[[change]]
        refit <- refits_list[[change_list$refit_filename]]
        if (!is.null(refit)) {
            
            param_infl <- model_grads$param_infls[[param]]
            
            param_infl$target_parameter
            target_index <- which(param_infl$target_parameter == model_grads$parameter_names)
            stopifnot(is.numeric(target_index))

            # beta_orig <- model_grads$betahat[target_index]
            # se_orig <- model_grads$se[target_index]
            beta_orig <- model_grads$model_fit$param[target_index]
            se_orig <- model_grads$model_fit$se[target_index]

            beta_rerun <- refit$betahat[target_index]
            se_rerun <- refit$se[target_index]
            
            
            beta_pred <- beta_orig + PredictChange(param_infl$param, change_list$amis)
            se_pred <- se_orig + PredictChange(param_infl$se, change_list$amis)
            
            desc_df <- data.frame(
                param=param,
                description=change,
                num_removed=change_list$n_drop,
                prop_removed=change_list$prop_drop)
            
            refits_df <- bind_rows(
                refits_df,
                bind_rows(
                    data.frame(beta=beta_orig, se=se_orig, fit="orig"),
                    data.frame(beta=beta_rerun, se=se_rerun, fit="rerun"),
                    data.frame(beta=beta_pred, se=se_pred, fit="pred")) %>%
                    bind_cols(desc_df)
            )
        }
    }
}


# refits_wide_df <-
#     refits_df %>% 
#     pivot_wider(id_cols=-c(beta, se, fit), 
#                 names_from=fit, values_from=c(beta, se))

refits_wide_df <-
  refits_df %>% 
  pivot_wider(names_from=fit, values_from=c(beta, se))

refit_filename <- file.path(
    base_dir, "data/", 
    sprintf("microcredit_vb_%s_refit_results.Rdata", analysis_desc))

print(sprintf("Saving to %s", refit_filename))
save(refits_df, refits_wide_df, amis_list, file=refit_filename)



#############################################################
# Exploratory and sanity checking graphs


if (FALSE) {
    warnings(1)
    ggplot(refits_wide_df) +
        geom_point(aes(x=beta_pred, y=beta_rerun, color=num_removed), size=2) +
        facet_grid( ~ param) +
        geom_abline(aes(slope=1, intercept=0))

    filter(refits_wide_df, str_detect(param, "log_sd_tau")) %>%
        mutate(direction=str_sub(description, 1, 3)) %>%
        ggplot() +
        geom_line(aes(x=num_removed, group=direction, color="Predicted", y=beta_pred - beta_orig), lwd=1.2) +
        geom_line(aes(x=num_removed, group=direction, color="Actual (rerun)", y=beta_rerun - beta_orig), lwd=1.2) +
        geom_point(aes(x=num_removed, group=direction, color="Predicted", y=beta_pred - beta_orig), size=2) +
        geom_point(aes(x=num_removed, group=direction, color="Actual (rerun)", y=beta_rerun - beta_orig), size=2) +
        facet_grid(param ~ direction) +
        ggtitle("Changes in log_sd_tau[*] after removing N most influential points.") +
        xlab("Number of points dropped") + ylab("Changes in log_sd_tau")

}




