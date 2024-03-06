#!/usr/bin/env Rscript

################################################################################
# Process the output of the notebooks
# microcredit_vb.ipynb and sensitivity.ipynb, and the R script run_refits.R
# into a format that we can use in the paper.

library(optparse)

if (interactive()) {
  git_repo_loc <- system("git rev-parse --show-toplevel", intern=TRUE)
  num_advi_draws <- 30
  python_path <- "/home/rgiordan/miniforge3/envs/amip-paper-2023/bin/python"
} else {
  option_list <- list(
    make_option("--git_repo_loc", type="character", default=NULL, 
                help="Paper github repo", metavar="character"),
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


library(ggplot2)
library(tidyverse)
library(zaminfluence)


# Set this to your path.
#git_repo_dir <- "/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench/"
git_repo_dir <- git_repo_loc
analysis_dir <- file.path(git_repo_dir, "examples/microcredit_mixture/")
paper_dir <- file.path(git_repo_dir, "writing/output/applications_data/microcredit_mixture/")

output_filename <- file.path(paper_dir, "microcredit_mixture_results.Rdata")

# Set these to the values for your particular analysis.
use_smuber <- "False"
#num_advi_draws <- opt$num_draws
analysis_desc <- sprintf("draws%d_smuber%s", num_advi_draws, use_smuber)

initial_fit_filename <- file.path(
    analysis_dir, "data/", 
    sprintf("microcredit_vb_%s_combined_results.Rdata", analysis_desc))

refit_filename <- file.path(
    analysis_dir, "data/", 
    sprintf("microcredit_vb_%s_refit_results.Rdata", analysis_desc))

load(initial_fit_filename)
load(refit_filename)


# Create a table with the same formatting as produced by zaminfluence
study_case_df <- data.frame()
for (s in c(0, 1)) {
    direction <- if(s == 1) "+" else "-"
    study_case_df <- bind_rows(
        study_case_df,
        data.frame(
            param_name=sprintf("tau[%d]", s),
            study_case=sprintf("$\\tau_{%s}$",
                               direction), stringsAsFactors=FALSE),
        data.frame(
            param_name=sprintf("log_sd_tau[%d]", s),
            study_case=sprintf("$\\log \\sigma_{\\tau_{%s}}$",
                               direction), stringsAsFactors=FALSE))
}


# names(reruns_zam_df)
# head(reruns_zam_df)

sig_num_ses <- qnorm(0.975)
all_results_df <- 
    refits_df %>%
    rename(param_name=param) %>% 
    rename(param=beta,
           n_drop=num_removed,
           prop_drop=prop_removed,
           analysis=fit) %>% 
    mutate(param_mzse=param - sig_num_ses * se,
           param_pzse=param + sig_num_ses * se,
           description=sub("sig$", "significance", description),
           description=sub("both", "sign and significance", description),
           study_case=ordered(
               param_name,
               levels=study_case_df$param_name,
               labels=study_case_df$study_case)) %>%
    pivot_longer(cols=c(param, param_mzse, param_pzse, se),
                 names_to="metric", values_to="value")

base_df <- filter(all_results_df, analysis=="orig")
results_df <- filter(all_results_df, analysis!="orig")

# names(reruns_zam_df)
# head(reruns_zam_df)
# names(results_df)

advi_time <- save_list[["advi_time"]]
#advi_mcmc_comparison_df <- save_list[["comp_df"]]

####################

print(sprintf("Saving to %s", output_filename))
save(base_df, results_df, num_advi_draws, advi_time,
     file=output_filename)

print("Success!")

# save(base_df, results_df,
#      num_advi_draws, advi_time, advi_mcmc_comparison_df,
#      file=output_filename)

# save(refits_df, refits_wide_df, results_df,
#      num_advi_draws, advi_mcmc_comparison_df, 
#      advi_time,
#      file=output_filename)