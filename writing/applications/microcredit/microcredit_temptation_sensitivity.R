#!/usr/bin/env Rscript
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

library(zaminfluence)
library(tidyverse)

base_dir <- git_repo_loc

# you should change this too if you want to write to a different location
output_dir <- file.path(base_dir, "writing/output/applications_data/microcredit")


################
# Load the data

load(file.path(output_dir, "microcredit_temptation_data.Rdata"))

# To prevent confusion epecially when rewriting later in the loop, give it a better name
full_Temptation_data <- data
rm(data)
N <- dim(full_Temptation_data)[1]
K <- length(unique(full_Temptation_data$site))

# Name the regression model.
reg_form <- formula("temptation ~ treatment + 1")

# create storage
results_list <- list()

# assess the sensitivity!
for(i in 1:K){
  data <- filter(full_Temptation_data, site == i) %>% mutate(row=1:n())
  reg_fit <- lm(data = data, formula = reg_form, x=TRUE, y=TRUE)

  # Get influence.
  model_grads <-
    ComputeModelInfluence(reg_fit) %>%
    AppendTargetRegressorInfluence("treatment")
  
  signals <- GetInferenceSignals(model_grads)
  reruns <- RerunForSignals(signals, model_grads)
  preds <- PredictForSignals(signals, model_grads)
  
  rerun_df <- 
    GetSignalsAndRerunsDataframe(signals, reruns, model_grads) %>%
    mutate(site=study_country[i])
  preds_df <- 
    GetSignalsAndRerunsDataframe(signals, preds, model_grads) %>%
    mutate(site=study_country[i])
  base_df <- 
    GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls) %>%
    mutate(site=study_country[i])
  
  results_list[[i]] <- list(rerun_df=rerun_df, base_df=base_df, reg_fit=reg_fit)
}


results_df <- bind_rows(
  map_df(results_list, ~ .$rerun_df) %>% mutate(analysis="rerun"),
  map_df(results_list, ~ .$preds_df) %>% mutate(analysis="prediction"))

base_df <- map_df(results_list, ~ .$base_df) %>% mutate(analysis="base")

save(results_df, base_df, file=file.path(output_dir, "microcredit_temptation_results.Rdata"))

print("Success!")
