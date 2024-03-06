# Reading Angelucci and Di Giorgi 2009 data for table 1
# Rachael Meager
# Dec 2018 

### Notes

### Preliminaries 

# install and load packages
#setwd("/home/rgiordan/Documents/git_repos/VariationalBayesPythonWorkbench/Models/Rachael/")
base_dir <- "/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench"

installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM", "MASS","quantreg",
                  "boot", "foreign")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}

# set the project WD first on your own console 
# on mine it is setwd("/Users/rachaelmeager/Dropbox/research work/research-work-LSE/automatic-robustness/")

# now read data 

#data <- read.dta("econ-examples/angelucci-digiorgi-data-code/table1.dta")
data <- read.dta(file.path(base_dir, "./examples/cash_transfers/angelucci-digiorgi-data-code/table1.dta"))

# some exporatory experiments

sum(is.na(data$treatp * data$treatnp)) # ok good



# now try to recreate table 1 for p, but just not clustering the SEs yet 
lm(Cind ~ treatp, data = data[data$t==10,])
lm(Cind ~ treatp, data = data[data$t==9,])
lm(Cind ~ treatp, data = data[data$t==8,])

# try to recreate table 1 for np, but just not clustering the SEs yet 
lm(Cind ~ treatnp, data = data[data$t==10,])
lm(Cind ~ treatnp, data = data[data$t==9,])
lm(Cind ~ treatnp, data = data[data$t==8,])

# not quite the same because didn't truncate, let's try

lm(Cind ~ treatp, data = data[data$t==10 & data$Cind <10000,])
lm(Cind ~ treatp, data = data[data$t==9 & data$Cind <10000,])
lm(Cind ~ treatp, data = data[data$t==8 & data$Cind <10000,])
lm(Cind ~ treatnp, data = data[data$t==10 & data$Cind <10000,])
lm(Cind ~ treatnp, data = data[data$t==9 & data$Cind <10000,])
lm(Cind ~ treatnp, data = data[data$t==8 & data$Cind <10000,])

# we've got it! 

# now let us bring in the controls

lm_treatnp_controls <- ( Cind ~ treatnp + avg_shock + hhhage + factor(hhhsex) + factor(p16) + factor(hhhalpha) + factor(hhhspouse) + yycali_1 + factor(region_1) + indice + hectareas + vhhnum )
lm_treatp_controls <- ( Cind ~ treatp + avg_shock + hhhage + factor(hhhsex) + factor(p16) + factor(hhhalpha) + factor(hhhspouse) + yycali_1 + factor(region_1) + indice + hectareas + vhhnum )


lm(lm_treatnp_controls,data = data[data$t==10 & data$Cind <10000,] )
lm(lm_treatnp_controls,data = data[data$t==9 & data$Cind <10000,] )
lm(lm_treatnp_controls,data = data[data$t==8 & data$Cind <10000,] )

lm(lm_treatp_controls,data = data[data$t==10 & data$Cind <10000,] )
lm(lm_treatp_controls,data = data[data$t==9 & data$Cind <10000,] )
lm(lm_treatp_controls,data = data[data$t==8 & data$Cind <10000,] )

# there are still some small discrepancies in the point estimates here, not sure why

# anyway, the only thing I havent done here is to cluster the standard errors and use the robust SE estimator 

dim(data)
dim(data[data$Cind <10000,])

saveRDS(data, file.path(base_dir, "examples/cash_transfers/processed_CT_data.rds"))
