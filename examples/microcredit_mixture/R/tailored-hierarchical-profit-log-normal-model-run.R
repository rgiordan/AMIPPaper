#!/usr/bin/env Rscript
# Command:
# 

# microcredit tailored hierarchical pdf for profit with lognormal tails
# Rachael Meager
# April 2016 

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

### Notes

# This file runs Stan on the microcredit model.

### Preliminaries 

repo_dir <- git_repo_loc
base_dir <- file.path(repo_dir, "/examples/microcredit_mixture/")
setwd(base_dir)


# 
# # install and load packages
# installPackage <- function(p) {
#     if (!is.element(p, installed.packages()[,1])) {
#         print(paste('installing', p))
#         install.packages(p, dep = TRUE)
#     } else {
#         print(paste(p, 'is already installed'))
#     }
#     require(p, character.only = TRUE) # test if it loads
# }

# installation_needed  <- FALSE
# loading_needed <- TRUE
# package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
#                   "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM")
# if(installation_needed){for (p in packagelist) {installPackage(p, repos='http://cran.us.r-project.org')}}
# if(loading_needed){lapply(package_list, require, character.only = TRUE)}

library(tidyverse)
library(rstan)
library(xtable)

# configure rstan

rstan_options(auto_write = TRUE)

#cores <- as.numeric(Sys.getenv('NSLOTS'))
cores <- as.numeric(6)

options(mc.cores = cores)

# load data
data <- read.csv(file=file.path(
  repo_dir, "examples/microcredit_mixture/R/data/microcredit_project_data_cleaned.csv"))


#stan needs these declared as data
N <- length(data$profit) # number of draws from the distribution / dimensionality of data vector
M <- 3 # mixture components
K <- length(unique(data$site)) # sites
P <- 2 # dimensions of X 
X <- cbind(rep(1,N), data$treatment)


cat <- rep(NA,N) # storage
for(i in 1:N){
  if(data$profit[i] < 0){ cat[i] <- 1  }
  else if(identical(data$profit[i],0)){cat[i] <- 2}
  else{cat[i] <- 3}
}


data_split <- list( data[cat==1,],data[cat==3,])

# now STAN MODEL! 


fit_data <- list(M = 3,
                 N = N,
                 K = K,
                 P = 2, 
                 Q = 10,
                 cat = cat,
                 site = data$site,
                 x = X, # This is used for the category prediction and is redudant.   
                 treatment_neg = data_split[[1]]$treatment,
                 y_neg = -data_split[[1]]$profit,
                 treatment_pos = data_split[[2]]$treatment,
                 y_pos = data_split[[2]]$profit,
                 N_neg = length(data_split[[1]]$profit),
                 N_pos = length(data_split[[2]]$profit),
                 site_neg = data_split[[1]]$site,
                 site_pos = data_split[[2]]$site,
                 quantile_vec = seq(0.05, 0.95, .1))


mcmc_time <- Sys.time()
stan_fit <- stan(file.path(base_dir, "R/stan-code/tailored-hierarchical-pdf-log-normal.stan"), 
                 iter = 4000, chains = 6, data = fit_data,
                 control = list(adapt_delta = 0.99, max_treedepth = 15))
mcmc_time <- Sys.time() - mcmc_time

print(stan_fit)

stan_fit_summary <- summary(stan_fit)
stan_fit_table <- xtable(stan_fit_summary$summary)

saveRDS(stan_fit_summary,
        file = "tailored_hierarchical_pdf_microcredit_output_4000_iters.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

save.image("output/microcredit_profit_lognormal_tailored_hierarchical_pdf_output_4000_iters.RData")

#sink("output/microcredit_profit_lognormal_tailored_hierarchical_output_table_4000.txt")
#print(stan_fit)
#dev.off()

print("Success!")






