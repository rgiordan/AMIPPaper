## Initial Script to Replicate Tables

library(dplyr)
library(broom)
library(sandwich)
library(purrr)
library(readr)
library(AER)

# This path is used in the mutliple testing github repo - need to change for other use.
setwd("/home/rgiordan/Documents/git_repos/VariationalBayesPythonWorkbench/Models/Rachael/oregon")

df <- read_csv("Data/OHI/OHIE_Public_Use_Files/OHIE_Data/data_for_analysis_R.csv")


##### Functions ####


# These functions just provide a wrapper to lm or ivreg with a few additions
# such as clustering on household and weighting the sample since the tables
# all use essentially the same specification. Rather than output an lm summary
# we use broom::tidy() to produce tidy regression output.


##### ITT
replicate_QJE_table_ITT <- function(outcome, dataset){
  dataset <- dataset %>% 
    filter(sample_12m_resp == 1) # This conditional seemed to be everywhere in the do file.
  
  # Creating model formula
  flma <- as.formula(paste0(outcome, "~ . + 0  - household_id"))
  # dplyr NSE
  # Only controls, treatment and outcome
  regression_df <- dataset %>%
    select_at(vars(treatment, 
                   household_id, # Need these for clustering SEs
                   starts_with("ddd"), # These are the fixed effects
                   outcome)) 
  # Hardcoded weight according to Stata do file for survey data.
  # RJG: weight_12m is "Final 12-month weights" from raw_dta_surveys[[3]]
  model <- lm(data = regression_df,
              formula = flma,
              weights = dataset$weight_12m)
  # Cluster on household
  clustered_SE <- vcovCL(model,
                         cluster = ~household_id) %>% 
    diag() %>% 
    sqrt()
  # Adding clustered standard errors, unadjusted p values and tidying
  tidy_model <- tidy(model,
                     quick = TRUE) %>% # quick means we compute our own SEs
    mutate(se = clustered_SE,
           tstat = estimate / se,
           pval = 2*pt(-abs(tstat), df = model$df))
  return(tidy_model)
}

##### LATE
replicate_QJE_table_LATE <- function(outcome, dataset){
  dataset <- dataset %>% 
    filter(sample_12m_resp == 1)
  # Creating model formula - removing instrument from reduced form
  flma <- as.formula(paste0(
    outcome,
    "~ . + 0 - household_id - treatment | . + 0 - household_id - ohp_all_ever_survey"))
  # This syntax is a bit funky. We use lm y ~ . to pick up all the FEs and then just
  # manually ensure the household_id isn't included (we need it for clustering later)
  # and the instrument appears in the right place.
  
  
  # Only controls, treatment and outcome
  regression_df <- dataset %>% 
    select_at(vars(treatment, # As above in ITT
                   household_id,
                   starts_with("ddd"),
                   outcome,
                   ohp_all_ever_survey)) # Instrument for receiving the OHP offer or something
  
  # Hardcoded weight according to Stata do file for survey data.
  model <- ivreg(data = regression_df,
                 formula = flma,
                 weights = dataset$weight_12m)
  
  # Cluster on household
  clustered_SE <- vcovCL(model, # This is a pet hate. ivreg won't play nicely with ~household_id
                         cluster = regression_df$household_id) %>% # specify it this way instead
    diag() %>% 
    sqrt()
  
  # Adding clustered standard errors, unadjusted p values and tidying
  tidy_model <- tidy(model) %>%
    select(term,
           estimate) %>% 
    mutate(se = clustered_SE,
           tstat = estimate / se,
           pval = 2*pt(-abs(tstat), df = model$df))
  return(tidy_model)
}


#### Table 5- any ####


# RJG: Why 12m?
ITT_table_5_any <- c("rx_any_12m",
                 "doc_any_12m",
                 "er_any_12m",
                 "hosp_any_12m") %>% 
  map_df(~(replicate_QJE_table_ITT(.x,
                             dataset = df) %>% 
        mutate(outcome = .x))) %>% 
  filter(term == "treatment")


LATE_table_5_any <-  c("rx_any_12m",
                   "doc_any_12m",
                   "er_any_12m",
                   "hosp_any_12m") %>% 
  map_df(~(replicate_QJE_table_LATE(.x,
                                     dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "ohp_all_ever_survey")

ITT_table_5_any %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), signif, 2) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "ITT Table V (any)")


LATE_table_5_any %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), round, 3) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "LATE Table V (any)")



##### Table 5 - number #####

ITT_table_5_number <- c("rx_num_mod_12m",
                        "doc_num_mod_12m",
                        "er_num_mod_12m",
                        "hosp_num_mod_12m") %>% 
  map_df(~(replicate_QJE_table_ITT(.x,
                                   dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "treatment")

LATE_table_5_number <- c("rx_num_mod_12m",
                         "doc_num_mod_12m",
                         "er_num_mod_12m",
                         "hosp_num_mod_12m") %>% 
  map_df(~(replicate_QJE_table_LATE(.x,
                                    dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "ohp_all_ever_survey")




ITT_table_5_number %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), signif, 2) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "ITT Table V (number)")

LATE_table_5_number %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), round, 3) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "LATE Table V (number)")




##### Table 6  - Compliance with Preventative Care ####

ITT_table_6 <- c("chl_chk_bin_12m",
                 "dia_chk_bin_12m",
                 "mam_chk_bin_12m",
                 "pap_chk_bin_12m") %>% 
  map_df(~(replicate_QJE_table_ITT(.x,
                                   dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "treatment")

LATE_table_6 <- c("chl_chk_bin_12m",
                  "dia_chk_bin_12m",
                  "mam_chk_bin_12m",
                  "pap_chk_bin_12m") %>% 
  map_df(~(replicate_QJE_table_LATE(.x,
                                    dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "ohp_all_ever_survey")


ITT_table_6 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), signif, 2) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "ITT Table VI")

LATE_table_6 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), round, 3) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "LATE Table VI")



##### Table 8 - Finance ####

ITT_table_8 <- c("cost_any_oop_12m",
                 "cost_any_owe_12m",
                 "cost_borrow_12m",
                 "cost_refused_12m") %>% 
  map_df(~(replicate_QJE_table_ITT(.x,
                                     dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "treatment")


LATE_table_8 <- c("cost_any_oop_12m",
                  "cost_any_owe_12m",
                  "cost_borrow_12m",
                  "cost_refused_12m") %>% 
  map_df(~(replicate_QJE_table_LATE(.x,
                                      dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "ohp_all_ever_survey")




ITT_table_8 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), signif, 2) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "ITT Table VIII")

LATE_table_8 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), round, 3) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "LATE Table VIII")
  



##### Table 9 - Health #####

ITT_table_9 <- c("health_genflip_bin_12m",
                 "health_notpoor_12m",
                 "health_chgflip_bin_12m",
                 "notbaddays_tot_12m",
                 "notbaddays_phys_12m",
                 "notbaddays_ment_12m",
                 "nodep_screen_12m") %>% 
  map_df(~(replicate_QJE_table_ITT(.x,
                                     dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "treatment")

LATE_table_9 <- c("health_genflip_bin_12m",
                  "health_notpoor_12m",
                  "health_chgflip_bin_12m",
                  "notbaddays_tot_12m",
                  "notbaddays_phys_12m",
                  "notbaddays_ment_12m",
                  "nodep_screen_12m") %>% 
  map_df(~(replicate_QJE_table_LATE(.x,
                                      dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "ohp_all_ever_survey")
  

ITT_table_9 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), signif, 2) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "ITT Table IX")

LATE_table_9 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), round, 3) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "LATE Table IX")



##### Table 10 ####
ITT_table_10 <- c("usual_clinic_12m",
                  "usual_doc_12m",
                  "needmet_med_12m",
                  "needmet_rx_12m",
                  "not_er_noner_12m") %>% 
  map_df(~(replicate_QJE_table_ITT(.x,
                                   dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "treatment")

LATE_table_10 <- c("usual_clinic_12m",
                   "usual_doc_12m",
                   "needmet_med_12m",
                   "needmet_rx_12m",
                   "not_er_noner_12m") %>% 
  map_df(~(replicate_QJE_table_LATE(.x,
                                    dataset = df) %>% 
             mutate(outcome = .x))) %>% 
  filter(term == "ohp_all_ever_survey")

ITT_table_10 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), signif, 2) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "ITT Table X")

LATE_table_10 %>% 
  select(-term) %>% 
  mutate_at(vars(estimate,
                 se,
                 tstat), round, 3) %>% 
  gt::gt() %>% 
  gt::fmt_scientific(columns = "pval") %>% 
  gt::tab_header(title = "LATE Table X")




##### Quick rstanarm comparison ####
# 
# 
# library(rstanarm)
# options(mc.cores = 4)
# ITT_table_5_num_rx_bayes <- stan_glmer(rx_num_mod_12m ~ treatment | numhh_list,
#                                      data = df %>% 
#                                        sample_frac(0.25)) 
# # Need to use 25% of sample to fit in decent time on my laptop.
# 
# library(tidybayes)
# model_draws <- ITT_table_5_num_rx_bayes %>%
#   spread_draws(b[term, numhh_list])
#   
# 
# 
# model_draws %>% 
#   median_qi() %>% 
#   mutate(frequentist = ifelse(term == "treatment",
#                               ITT_table_5_number %>%
#                                 filter(outcome == "rx_num_mod_12m") %>%
#                                 select(estimate) %>%
#                                 pull(),
#                               NA)) %>% 
#   ggplot(aes(x = b,
#              xmin = .lower,
#              xmax = .upper,
#              y = numhh_list,
#              colour = numhh_list)) +
#   ggstance::geom_pointrangeh() +
#   facet_wrap(~term, scales = "free") +
#   theme_minimal() +
#   geom_vline(aes(xintercept = frequentist), linetype = "longdash") +
#   labs(caption = "Dashed line indicates frequentist pooled estimate. \n
#        N.B. Ignoring survey fixed effects.",
#        title = "No. prescription drugs after treatment",
#        subtitle = "Using 25% of sample and ignoring survey FEs")






