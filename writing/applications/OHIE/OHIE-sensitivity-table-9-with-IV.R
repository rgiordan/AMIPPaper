#!/usr/bin/env Rscript

library(sandwich)
library(AER)
library(zaminfluence)
library(readr)
library(tidyverse)


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


# It should be enough to change this variable to point to your own installation.
#base_dir <- "/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench"
base_dir <- git_repo_loc

# you should change this too if you want to write output
output_dir <- file.path(base_dir, "writing/output/applications_data/ohie")
data_dir <- file.path(base_dir, "examples/oregon/Data/OHI/OHIE_Public_Use_Files/OHIE_Data")

################
# Load the data
# this was pre-processed in the prepare_OHIE_data.R script

orig_df <- read_csv(file.path(data_dir, "data_for_analysis_R.csv"))

table_9_outcomes <-
  c("health_genflip_bin_12m",
    "health_notpoor_12m",
    "health_chgflip_bin_12m",
    "notbaddays_tot_12m",
    "notbaddays_phys_12m",
    "notbaddays_ment_12m",
    "nodep_screen_12m")

# Loop over outcomes.
reg_results_list <- list()
i <- 1
for(i in 1:length(table_9_outcomes)) {
  print(sprintf("Running for outcome %d", i))
  outcome <- table_9_outcomes[i]

  df <- orig_df %>%
    filter(sample_12m_resp == 1) %>%
    mutate(row=1:n())

  # Creating model formula
  flma <- as.formula(paste0(outcome, "~ . + 0  - household_id"))
  reg_formula <- flma
  target_regressor <- "treatment"

    # dplyr NSE
  # Only controls, treatment and outcome
  regression_df <- df %>%
    select_at(vars(treatment,
                   household_id, # Need these for clustering SEs
                   starts_with("ddd"), # These are the fixed effects
                   !!outcome))

  # Note: weight_12m is "Final 12-month weights" from raw_dta_surveys[[3]]
  # to get the right number of rows we have to staple and un-staple the weights, i believe
  regression_df$weights <- df$weight_12m
  regression_df <- regression_df[complete.cases(regression_df), ]
  complete_cases_weight_12m <- regression_df$weights
  regression_df$weights <- NULL # this is to make the formula with "." work

  # re-index the cluster index
  regression_df$household_id <- regression_df$household_id - min(regression_df$household_id)

  # Hardcoded weight according to Stata do file for survey data.
  model <- lm(data = regression_df, weights = complete_cases_weight_12m,
              formula = reg_formula, x = TRUE, y = TRUE)

  summary(model)

  cat("Regression influence...")
  model_grads <-
    ComputeModelInfluence(model, se_group=regression_df$household_id) %>%
    AppendTargetRegressorInfluence("treatment")

  signals <- GetInferenceSignals(model_grads)
  reruns <- RerunForSignals(signals, model_grads)
  preds <- PredictForSignals(signals, model_grads)

  rerun_df <-
    GetSignalsAndRerunsDataframe(signals, reruns, model_grads) %>%
    mutate(outcome=outcome, method="regression")
  preds_df <-
    GetSignalsAndRerunsDataframe(signals, preds, model_grads) %>%
    mutate(outcome=outcome, method="regression")
  base_df <-
    GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls) %>%
    mutate(outcome=outcome, method="regression")

  results <- list(
    rerun_df=rerun_df,
    preds_df=preds_df,
    base_df=base_df,
    lm_result=model
  )
  reg_results_list[[i]] <- results
}


iv_results_list <- list()
for(i in 1:length(table_9_outcomes)) {
  print(sprintf("Running for outcome %d", i))
  outcome <- table_9_outcomes[i]

  ## can we do some IV
  df <- orig_df %>%
    filter(sample_12m_resp == 1) %>%
    mutate(row=1:n())

  # Creating model formula - removing instrument from reduced form
  flma <- as.formula(paste0(
    outcome,
    "~ . + 0 - household_id - treatment | . + 0 - household_id - ohp_all_ever_survey"))
  # This syntax is a bit funky. We use lm y ~ . to pick up all the FEs and then just
  # manually ensure the household_id isn't included (we need it for clustering later)
  # and the instrument appears in the right place.

  # Only controls, treatment, instrument and outcome
  regression_df <- df %>%
    select_at(vars(treatment, # As above in ITT
                   household_id,
                   starts_with("ddd"),
                   !!outcome,
                   ohp_all_ever_survey)) # Instrument for receiving the OHP offer or something

  # Hardcoded weight according to Stata do file for survey data.
  regression_df$weights <- df$weight_12m
  regression_df <- regression_df[complete.cases(regression_df),]
  complete_cases_weight_12m <- regression_df$weights
  regression_df$weights <- NULL
  # Hardcoded weight according to Stata do file for survey data.
  model <- ivreg(data = regression_df,weights = complete_cases_weight_12m,
                 formula = flma, y = TRUE, x = TRUE)

  # now we must re-base our household ID
  regression_df$household_id <- regression_df$household_id - min(regression_df$household_id)

  model_grads <-
    ComputeModelInfluence(model, se_group=regression_df$household_id) %>%
    AppendTargetRegressorInfluence("ohp_all_ever_survey")

  signals <- GetInferenceSignals(model_grads)
  reruns <- RerunForSignals(signals, model_grads)
  preds <- PredictForSignals(signals, model_grads)

  rerun_df <-
    GetSignalsAndRerunsDataframe(signals, reruns, model_grads) %>%
    mutate(outcome=outcome, method="iv")
  preds_df <-
    GetSignalsAndRerunsDataframe(signals, preds, model_grads) %>%
    mutate(outcome=outcome, method="iv")
  base_df <-
    GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls) %>%
    mutate(outcome=outcome, method="iv")

  cat("done.\n")

  results <- list(
    rerun_df=rerun_df,
    preds_df=preds_df,
    base_df=base_df,
    ivreg_result=model
  )
  iv_results_list[[i]] <- results
}

results_df <- bind_rows(
  map_df(reg_results_list, ~ .$rerun_df) %>% mutate(analysis="rerun"),
  map_df(reg_results_list, ~ .$preds_df) %>% mutate(analysis="prediction"),
  map_df(iv_results_list, ~ .$rerun_df) %>% mutate(analysis="rerun"),
  map_df(iv_results_list, ~ .$preds_df) %>% mutate(analysis="prediction")
)

base_df <- bind_rows(
  map_df(reg_results_list, ~ .$base_df) %>% mutate(analysis="base"),
  map_df(iv_results_list, ~ .$base_df) %>% mutate(analysis="base")
)


output_full_path <- file.path(output_dir, "OHIE_results.Rdata")
print(sprintf("Saving to %s", output_full_path))

save(results_df, base_df, file=output_full_path)

print("Success!")
