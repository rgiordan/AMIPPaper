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

base_dir <- git_repo_loc

# you should change this too if you want to write
output_dir <- file.path(base_dir, "writing/output/applications_data/microcredit")

################
# Load the data

load(file.path(base_dir, "examples/microcredit/microcredit_project_data.RData"))
# IF for some reason the above does not work then sometimes this does 
#load("/Users/rachaelmeager/AdversarialInfluenceWorkbench/examples/microcredit/microcredit_project_data.RData")

if( exists("USD_convert_to_2009_dollars")!=TRUE){stop("OLD DATA FILE LOADED! RELOAD CORRECT FILE!")}
# This preps data so that the model can be written efficiently - it's not necessary but
# it's the best way I know to code it in STAN
site <- c(
    angelucci_indicator, 
    attanasio_indicator, 
    augsberg_indicator,
    banerjee_indicator,
    crepon_indicator, 
    karlan_indicator, 
    tarozzi_indicator)
profit <- c(
    angelucci_profit, 
    attanasio_profit, 
    augsberg_profit,
    banerjee_profit,
    crepon_profit, 
    karlan_profit, 
    tarozzi_profit)
treatment <- c(
    angelucci_treatment, 
    attanasio_treatment, 
    augsberg_treatment,
    banerjee_treatment, 
    crepon_treatment, 
    karlan_treatment, 
    tarozzi_treatment)

# Now we have to standardise any variables which are in local currency units to USD PPP
# per fortnight in 2009 dollars

expanded_profit_standardiser_USD_PPP_per_fortnight <-
    c( rep(the_profit_standardiser_USD_PPP_per_fortnight[1], length(angelucci_indicator)),
       rep(the_profit_standardiser_USD_PPP_per_fortnight[2], length(attanasio_indicator)),
       rep(the_profit_standardiser_USD_PPP_per_fortnight[3], length(augsberg_indicator)),
       rep(the_profit_standardiser_USD_PPP_per_fortnight[4], length(banerjee_indicator)),
       rep(the_profit_standardiser_USD_PPP_per_fortnight[5], length(crepon_indicator)),
       rep(the_profit_standardiser_USD_PPP_per_fortnight[6], length(karlan_indicator)),
       rep(the_profit_standardiser_USD_PPP_per_fortnight[7], length(tarozzi_indicator)))

#profit <- profit * expanded_profit_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(
    site=site,
    profit=profit * expanded_profit_standardiser_USD_PPP_per_fortnight,
    treatment=treatment)

# We gotta remove the NA values
data <- data[complete.cases(data),]

print(sprintf("Saving to %s", file.path(output_dir, "microcredit_profit_data.Rdata")))
save(data, study_country, file=file.path(output_dir, "microcredit_profit_data.Rdata"))


######################################################################
# The temptation analysis omits a couple of sites

# Now we have to standardise any variables which are in local currency units to USD PPP
# per fortnight in 2009 dollars

site <- c(
    angelucci_indicator, 
    attanasio_indicator, 
    augsberg_indicator,
    banerjee_indicator,
    crepon_indicator)
treatment <- c(
    angelucci_treatment, 
    attanasio_treatment, 
    augsberg_treatment,
    banerjee_treatment, 
    crepon_treatment)
temptation <- c(
    angelucci_temptation, 
    attanasio_temptation, 
    augsberg_temptation,banerjee_temptation,
    crepon_temptation)

expanded_temptation_standardiser_USD_PPP_per_fortnight <-
    c( rep(the_temptation_standardiser_USD_PPP_per_fortnight[1], length(angelucci_indicator)),
       rep(the_temptation_standardiser_USD_PPP_per_fortnight[2], length(attanasio_indicator)),
       rep(the_temptation_standardiser_USD_PPP_per_fortnight[3], length(augsberg_indicator)),
       rep(the_temptation_standardiser_USD_PPP_per_fortnight[4], length(banerjee_indicator)),
       rep(the_temptation_standardiser_USD_PPP_per_fortnight[5], length(crepon_indicator)))

#temptation <- temptation*expanded_temptation_standardiser_USD_PPP_per_fortnight

# bind everything into a data frame
data <- data.frame(
    site=site, 
    temptation=temptation*expanded_temptation_standardiser_USD_PPP_per_fortnight,
    treatment=treatment)

# We gotta remove the NA values
data <- data[complete.cases(data),]

print(sprintf("Saving to %s", file.path(output_dir, "microcredit_temptation_data.Rdata")))
save(data, study_country, file=file.path(output_dir, "microcredit_temptation_data.Rdata"))

print("Success!")