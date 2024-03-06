#!/usr/bin/env Rscript
# Postprocess the cash transfers data for the paper.

# This file consumes the output of read-cash-data-into-R.R and produces
# the files necessary to produce the graphs and tables in the paper.

library(optparse)

if (interactive()) {
  git_repo_loc <- system("git rev-parse --show-toplevel", intern=TRUE)
} else {
  option_list <- list(
    make_option("--git_repo_loc", type="character", default=NULL, 
                help="Paper github repo", metavar="character")
  )
  opt_parser <- OptionParser(option_list=option_list);
  opt <- parse_args(opt_parser);
  git_repo_loc <- opt$git_repo_loc
}
if (!dir.exists(git_repo_loc)) {
  stop(sprintf("Repository directory %s does not exist.", git_repo_loc))
}


library(tidyverse)
library(gridExtra)
library(sandwich)
library(AER)
library(zaminfluence)


# It should be enough to change this variable to point to your own installation.
base_dir <- git_repo_loc
#base_dir <- "/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench"
#base_dir <- "/Users/rachaelmeager/AdversarialInfluenceWorkbench"

# you should change this too if you want to write output
output_dir <- file.path(base_dir, "writing/output/applications_data/cash_transfers")


library(devtools)
# Optionally re-load zaminfluence while developing
#load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")


################
# Load the data

full_cash_data <- readRDS(
  file.path(base_dir, "examples/cash_transfers/processed_CT_data.rds"))


# Truncate as in the original paper,
# define a "treatment" variable according to the regression variable,
# and remove experiments where both treatnp and treatp were NA.
data <- full_cash_data %>%
  filter(Cind < 10000) %>%
  select(-treatd, -treatap) %>%
  filter(reg %in% c("treatnp", "treatp")) %>%
  mutate(treatment=case_when(
      reg == "treatnp" ~ treatnp,
      reg == "treatp" ~ treatp,
      TRUE ~ as.integer(NA)))

# Sanity check
stopifnot(!any(is.na(data$treatment)))

# generic reg formula, which will only work when the treatnp and treatp
# variables are renamed below
reg_form <- formula(
  Cind ~ treatment + avg_shock + hhhage + factor(hhhsex) + factor(p16) +
         factor(hhhalpha) + factor(hhhspouse) + yycali_1 + factor(region_1) +
         indice + hectareas + vhhnum )

# assess the sensitivity!
# Break it up into the two types of cases
# Define the order of the cases here
study_cases <- c()
for (reg in c("treatnp", "treatp")) {
  for (time in c(10, 9, 8)) {
    study_cases <- c(study_cases, sprintf("%s, t=%d", reg, time))
  }
}
stopifnot(setequal(study_cases, unique(data$case)))

# Run sensitivity for all the cases.
result_list <- list()
for(i in 1:length(study_cases)){
  study_case <- study_cases[i]
  print(sprintf("Running for %s", study_case))
  case_data <- filter(data, case == study_case) %>% select(-treatnp, -treatp)
  case_data <- case_data[complete.cases(case_data), ]
  
  reg_fit <- lm(data=case_data, formula=reg_form, x=TRUE, y=TRUE)

  RerunFun <- function(w) {
    # Use a custom re-running function in zaminfluence, since
    # some of the weight changes drop a level from a factor,
    # resulting in a (removable) singularity in the space of weights.
    #
    # When you pass the modified dataframe into R's lm function,
    # the missing level's indicator is automatically dropped.
    
    keep_rows <- abs(w) > 1e-8
    data_for_reg_new <- case_data[keep_rows, ]
    lm_rerun <- lm(data=data_for_reg_new, formula=reg_form)
    cluster_var <- vcovCL(lm_rerun,
                          cluster=data_for_reg_new$village,
                          type="HC0", cadjust=FALSE)
    return(ModelFit(fit_object=lm_rerun, 
                    num_obs=nrow(data_for_reg_new), 
                    param=lm_rerun$coefficients,
                    se=sqrt(diag(cluster_var)),
                    parameter_names=names(lm_rerun$coefficients), 
                    se_group=model_grads$model_fit$se_group[keep_rows]))
  }

  # Get influence and reruns.
  model_grads <-
    ComputeModelInfluence(reg_fit, se_group=case_data$village) %>%
    AppendTargetRegressorInfluence("treatment")
  
  signals <- GetInferenceSignals(model_grads)
  reruns <- RerunForSignals(signals, model_grads, RerunFun=RerunFun)
  preds <- PredictForSignals(signals, model_grads)

  rerun_df <- GetSignalsAndRerunsDataframe(signals, reruns, model_grads)
  preds_df <- GetSignalsAndRerunsDataframe(signals, preds, model_grads)
  base_df <- GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls)
  
  # Sanity check that the zaminfluence covariances match vcovCL. 
  vcov_se <- vcovCL(
    reg_fit,
    cluster=case_data$village, type="HC0", cadjust=FALSE) %>% diag() %>% sqrt()
  stopifnot(max(abs(model_grads$model_fit$se - vcov_se)) < 1e-8)

  results <- list()
  results$reg_fit <- reg_fit
  results$rerun_df <- rerun_df %>% mutate(study_case=study_case)
  results$preds_df <- preds_df %>% mutate(study_case=study_case)
  results$base_df <- base_df %>% mutate(study_case=study_case)
  result_list[[i]] <- results
}

results_df <- bind_rows(
  map_df(result_list, ~ .$rerun_df) %>% mutate(analysis="rerun"),
  map_df(result_list, ~ .$preds_df) %>% mutate(analysis="prediction"))

base_df <- map_df(result_list, ~ .$base_df) %>% mutate(analysis="base")

save(results_df, base_df, result_list,
     file=file.path(output_dir, "cash_transfers_results.Rdata"))

print("Success!")