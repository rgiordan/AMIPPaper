# Use this script to debug and edit the knit graphs without re-compiling in latex.

base_dir <- "/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench"
#base_dir <- "/Users/rachaelmeager/AdversarialInfluenceWorkbench"

paper_directory <- file.path(base_dir, "writing/output/")

knitr_debug <- FALSE # Set to true to see error output
cash_cache <- FALSE
ohie_cache <- FALSE
microcredit_cache <- FALSE

setwd(paper_directory)
source(file.path(paper_directory, "figures_knitr/initialize.R"))
source(file.path(paper_directory, "figures_knitr/load_data.R"))

# Now you can source individual files from the paper to see how they look
# without re-compiling.


knitr_ohie_file <- file.path(data_path, "ohie", "OHIE_results.Rdata")
knitr_ohie_file == output_full_path
cat(knitr_ohie_file)
cat(output_full_path)

ohie_env <- LoadIntoEnvironment(
  file.path(data_path, "ohie", "OHIE_results.Rdata"))


source("figures_knitr/OHIE/OHIE-table_results.R",
       echo=knitr_debug, print.eval=TRUE)

ohie_env$results_df %>% 
  filter(outcome == "health_genflip_bin_12m") %>%
  filter(metric == "param", method == "iv", target_qoi=="param") %>%
  select(value, target_param_name, n_drop, outcome, analysis)


# Running writing/applications/OHIE/OHIE-sensitivity-table-9-with-IV.R
results_df %>%
  filter(outcome == "health_genflip_bin_12m") %>%
  filter(metric == "param", method == "iv", target_qoi=="param") %>%
  select(value, target_param_name, n_drop, outcome, analysis)
