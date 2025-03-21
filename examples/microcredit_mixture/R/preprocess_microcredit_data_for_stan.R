#!/usr/bin/env Rscript

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

library(tidyverse)
#package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
#                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM")
#if(installation_needed){for (p in packagelist) {installPackage(p, repos='http://cran.us.r-project.org')}}
#if(loading_needed){lapply(package_list, require, character.only = TRUE)}

load(file.path(repo_dir, "examples/microcredit/microcredit_project_data.RData"))

if( exists("USD_convert_to_2009_dollars")!=TRUE){stop("OLD DATA FILE LOADED! RELOAD CORRECT FILE!")}
# This preps data so that the model can be written efficiently - it's not necessary but it's the best way I know to code it in STAN
site <- c(angelucci_indicator, attanasio_indicator, augsberg_indicator,banerjee_indicator, crepon_indicator, karlan_indicator, tarozzi_indicator)
profit <- c(angelucci_profit, attanasio_profit, augsberg_profit,banerjee_profit, crepon_profit, karlan_profit, tarozzi_profit)
treatment <- c(angelucci_treatment, attanasio_treatment, augsberg_treatment,banerjee_treatment, crepon_treatment, karlan_treatment, tarozzi_treatment)

# Now we have to standardise any variables which are in local currency units to USD PPP per fortnight in 2009 dollars


expanded_standardiser_USD_PPP_per_fortnight <- c( rep(the_profit_standardiser_USD_PPP_per_fortnight[1],length(angelucci_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[2],length(attanasio_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[3],length(augsberg_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[4],length(banerjee_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[5],length(crepon_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[6],length(karlan_indicator)),
                                                  rep(the_profit_standardiser_USD_PPP_per_fortnight[7],length(tarozzi_indicator)))

profit <- profit*expanded_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(site, profit, treatment)

# We gotta remove the NA values
data <- data[complete.cases(data),]

# create categorical allocations
N <- length(data$profit) # number of draws from the distribution / dimensionality of data vector
cat <- rep(NA,N) # storage
for(i in 1:N){
  if(data$profit[i] < 0){ cat[i] <- 1  }
  else if(identical(data$profit[i],0)){cat[i] <- 2}
  else{cat[i] <- 3}
}

data <- data.frame(data, cat)

# Save the data in csv format
write.csv(data, file=file.path(base_dir, "R/data/microcredit_project_data_cleaned.csv"), row.names=FALSE)
print("Success!")






