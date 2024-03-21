#!/usr/bin/env Rscript
# Reading Angelucci and Di Giorgi 2009 data for table 1
# Rachael Meager
# Dec 2018

### Notes

### Preliminaries

# set the project WD first on your own console
# It should be enough to change this variable to point to your own installation.

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

base_dir <- git_repo_loc

# now read data
library(foreign)
library(tidyverse)



filepath_data <- file.path(base_dir, "examples/cash_transfers/angelucci-digiorgi-data-code/table1.dta")
if (!file.exists(filepath_data)) {
  stop(sprintf("%s does not exist", filepath_data))
}
raw_data <- read.dta(filepath_data)

### CHECK THAT WE CAN REPLICATE THE PAPER

# some exporatory experiments

sum(is.na(raw_data$treatp * raw_data$treatnp)) # ok good


# now try to recreate table 1 for p, but just not clustering the SEs yet
# Note that if we do not truncat Cind then we don't match.
for (time in c(10, 9, 8)) {
    cat("\n------------------------\nt =", time, ":\n")
    print(lm(Cind ~ treatp, data = raw_data %>% filter(t == time, Cind < 10000)))
    print(lm(Cind ~ treatnp, data = raw_data %>% filter(t == time, Cind < 10000)))
}

# we've got it!

# now let us bring in the controls
controls <- paste(c(
    "avg_shock",
    "hhhage",
    "factor(hhhsex)",
    "factor(p16)",
    "factor(hhhalpha)",
    "factor(hhhspouse)",
    "yycali_1",
    "factor(region_1)",
    "indice",
    "hectareas",
    "vhhnum"), collapse=" + ")

lm_treatnp_controls <- formula(paste0("Cind ~ treatnp + ", controls))
lm_treatp_controls <- formula(paste0("Cind ~ treatp + ", controls))

for (time in c(10, 9, 8)) {
    cat("\n------------------------\nt =", time, ":\n")
    print(lm(lm_treatp_controls, data = raw_data %>% filter(t == time, Cind < 10000)))
    print(lm(lm_treatnp_controls, data = raw_data %>% filter(t == time, Cind < 10000)))
}

# there are still some small discrepancies in the point estimates here, not sure why

# now let us put it into cases which are easy to refer to in a loop later

# The data is such that only one of treatnp or treatp is not NA, so the following division of
# case_indicator is a partition.
stopifnot(sum((!is.na(raw_data$treatnp)) & (!is.na(raw_data$treatp))) == 0)

raw_data <- raw_data %>%
    mutate(reg=case_when(is.na(treatnp) & !is.na(treatp) ~ "treatp",
                         !is.na(treatnp) & is.na(treatp) ~ "treatnp",
                         TRUE ~ "error"),
           case=sprintf("%s, t=%d", reg, t))


processed_data_location <- file.path(base_dir, "examples/cash_transfers/processed_CT_data.rds")
saveRDS(raw_data, file=processed_data_location)

print("Success!")
