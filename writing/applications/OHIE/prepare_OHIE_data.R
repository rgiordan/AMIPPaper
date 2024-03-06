#!/usr/bin/env Rscript
# Replicating OHIE Data

# This script replicates prepare_data.do used in the QJE paper.
# The output is near idential - we use a different naming convention for the 
# fixed effects but use the same prefixes as the stata file so dplyr's
# starts_with("eee") will give the correct corresponding FEs for instance.
# Stata seems to expand FEs a little different to R so all our models are 
# saturated - setting intercept = 0 or dropping a set of FEs is required to run
# most regressions.

# ITT and LATE estimates from all tables of the paper replicate perfectly and I run
# a few basic unit tests throughout.

# The script returns the clean dataset as a tibble, alternatively you can save to
# a CSV.


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


write_to_disk <- TRUE

library(haven)
library(dplyr)
library(purrr)

data_dir <- file.path(
  git_repo_loc, "examples/oregon/Data/OHI/OHIE_Public_Use_Files/OHIE_Data")

## Loading DTAs
raw_dta_desc <- read_dta(file.path(data_dir, "oregonhie_descriptive_vars.dta"))
raw_dta_surveys <- c(
  file.path(data_dir, "oregonhie_survey0m_vars.dta"),
  file.path(data_dir, "oregonhie_survey6m_vars.dta"),
  file.path(data_dir, "oregonhie_survey12m_vars.dta")) %>% 
  map(read_dta)
raw_dta_state <- read_dta(file.path(data_dir, "oregonhie_stateprograms_vars.dta"))


## Merging to desc df
merge_df <- inner_join(raw_dta_desc,
                       raw_dta_surveys[[1]],
                       by = "person_id") %>% 
  inner_join(.,
             raw_dta_surveys[[2]],
             by = "person_id") %>% 
  inner_join(.,
             raw_dta_surveys[[3]],
             by = "person_id") %>% 
  inner_join(.,
             raw_dta_state,
             by = "person_id")

# Will throw an error if we've lost someone
testthat::test_that("matches ok", {
  testthat::expect_equal(nrow(raw_dta_desc), nrow(merge_df))
})


## Renaming Insurance Variables - funfunfunfufnunfufn
merge_df <- merge_df %>% 
  rename(
    ohp_all_ever_admin =  ohp_all_ever_matchn_30sep2009, 
    ohp_all_mo_admin =ohp_all_mo_matchn_30sep2009,   
    ohp_all_mo_survey =ohp_all_mo_firstn_30sep2009,   
    ohp_all_ever_survey =ohp_all_ever_firstn_30sep2009, 
    ohp_all_ever_survey0m=ohp_all_ever_firstn_survey0m,  
    ohp_all_ever_survey6m=ohp_all_ever_firstn_survey6m,  
    ohp_all_end_admin =ohp_all_end_30sep2009,         
    ohp_all_mo_survey0m =ohp_all_mo_firstn_survey0m,    
    ohp_all_mo_survey6m =ohp_all_mo_firstn_survey6m    
  ) %>% 
  mutate(
    ohp_all_end_survey = ohp_all_end_admin
  ) %>% 
  rename(
    ohp_std_ever_admin =ohp_std_ever_matchn_30sep2009, 
    # ohp_std_mo_admin =ohp_std_mo_matchn_30sep2009,   
    ohp_std_ever_survey =ohp_std_ever_firstn_30sep2009 
    # ohp_std_mo_survey0m=ohp_std_mo_firstn_survey0m,    
    # ohp_std_mo_survey6m=ohp_std_mo_firstn_survey6m   
  ) %>% 
  rename(
    postn_tanf_bin = tanf_ever_matchn_30sep2009,   
    prenany_tanf_bin = tanf_ever_prenotify07,        
    postn_survey12m_tanf_bin = tanf_ever_firstn_survey12m,   
    pren_survey12m_tanf_bin= tanf_ever_presurvey12m,       
    postn_survey12m_tanf_hh_amt = tanf_tot_hh_firstn_survey12m, 
    postn_tanf_hh_amt = tanf_tot_hh_30sep2009,        
    prenany_tanf_hh_amt = tanf_tot_hh_prenotify07,      
    pren_survey12m_tanf_hh_amt = tanf_tot_hh_presurvey12m     
  ) %>% 
  rename(
    postn_snap_bin = snap_ever_matchn_30sep2009,   
    pren_survey12m_snap_bin  = snap_ever_presurvey12m,       
    postn_snap_hh_amt = snap_tot_hh_30sep2009,        
    prenany_snap_hh_amt  = snap_tot_hh_prenotify07,      
    pren_survey12m_snap_hh_amt  = snap_tot_hh_presurvey12m,     
    # ohp_all_at_12m = ohp_all_at_survey12m,         
    prenany_snap_bin  = snap_ever_prenotify07,        
    postn_survey12m_snap_bin = snap_ever_firstn_survey12m,   
    postn_survey12m_snap_hh_amt = snap_tot_hh_firstn_survey12m 
  ) %>% 
  rename(
    zip_msa = zip_msa_list,
    draw = wave_survey0m,
    draw_survey_12m = wave_survey12m)

## Regression fixed effects

# 80% certain this imitates stata's xi command
options(na.action='na.pass') # NAs only in this bit
dummy_survey_draws <- merge_df %>% 
  mutate_at(vars(draw_survey_12m, numhh_list), as.factor) %>% 
  model.matrix(object = ~ 0 + draw_survey_12m:numhh_list) %>% 
  as_tibble() %>% 
  select(-`draw_survey_12m4:numhh_list3`,
         -`draw_survey_12m5:numhh_list3`,
         -`draw_survey_12m6:numhh_list3`,
         -`draw_survey_12m7:numhh_list3`)
# This mimics the prefix command in Stata although colname syntax isn't a perfect match
colnames(dummy_survey_draws) <- paste0("ddd", colnames(dummy_survey_draws))
options(na.action = "na.fail")


dummy_draw_lottery <- merge_df %>% 
  mutate(draw_lottery = factor(draw_lottery)) %>%
  select(draw_lottery) %>% 
  model.matrix(object = ~ 0 + draw_lottery) %>% 
  as_tibble()
colnames(dummy_draw_lottery) <- paste0("lll", colnames(dummy_draw_lottery))

dummy_draw_numhh <- merge_df %>% 
  select(numhh_list) %>%
  mutate(numhh_list = factor(numhh_list)) %>% 
  model.matrix(object = ~ 0 + numhh_list) %>% 
  as_tibble()
colnames(dummy_draw_numhh) <- paste0("nnn", colnames(dummy_draw_numhh))


dummy_df <- merge_df %>% 
  bind_cols(dummy_survey_draws,
            dummy_draw_lottery,
            dummy_draw_numhh)


# Expected number of rows?
testthat::test_that("correct rows", {
  testthat::expect_equal(
    list(merge_df,
         dummy_survey_draws,
         dummy_draw_lottery,
         dummy_draw_numhh,
         dummy_df) %>% 
      map(nrow) %>% 
      unlist(),
    rep(74922, 5)
  )
}) # YES!



# weights

dummy_df <- dummy_df %>% 
  mutate(noweight = 1,
         constant = 1,
         weight_0m = 1,
         weight_0m = ifelse(sample_0m != 1, NA, weight_0m),
         sample_0m_resp = (returned_0m == 1 & weight_0m != 0),
         sample_6m_resp = (returned_6m == 1 & weight_6m != 0))

# More FEs
options(na.action='na.pass')
dummy_draw_num_2 <- dummy_df %>% 
  select(draw, numhh_list) %>% 
  mutate_all(as.factor) %>% 
  model.matrix(object = ~ 0 + draw:numhh_list) %>% 
  as_tibble() %>% 
  select(-`draw4:numhh_list3`,
         -`draw5:numhh_list3`,
         -`draw6:numhh_list3`,
         -`draw7:numhh_list3`,
         -`draw8:numhh_list3`)
colnames(dummy_draw_num_2) <- paste0("eee", colnames(dummy_draw_num_2))
options(na.action = "na.fail")
dummy_df <- bind_cols(
  dummy_df,
  dummy_draw_num_2
)

# Other Variables




other_variable_df <- dummy_df %>% 
  mutate(health_poor_0m = (health_gen_0m == 1),
         health_poor_6m = (health_gen_6m == 1)) %>%
  mutate(
    health_genflip_bin_0m = case_when(
      health_gen_bin_0m == 1 ~ 0,
      health_gen_bin_0m == 0 ~ 1
  ),
  health_genflip_bin_6m = case_when(
    health_gen_bin_6m == 1 ~ 0,
    health_gen_bin_6m == 0 ~ 1
  ),
  health_notpoor_0m = case_when(
    health_poor_0m == 1 ~ 0,
    health_poor_0m == 0 ~ 1
  ),
  health_notpoor_6m = case_when(
    health_poor_6m == 1 ~ 0,
    health_poor_6m == 0 ~ 1
  ),
  health_chgflip_bin_0m = case_when(
    health_chg_bin_0m == 1 ~ 0,
    health_chg_bin_0m == 0 ~ 1
  ),
  health_chgflip_bin_6m = case_when(
    health_chg_bin_6m == 1 ~ 0,
    health_chg_bin_6m == 0 ~ 1
  ),
  notbaddays_tot_0m = 30 - baddays_tot_0m,
  notbaddays_phys_0m = 30 - baddays_phys_0m,
  notbaddays_ment_0m = 30 - baddays_ment_0m,
  notbaddays_tot_6m = 30 - baddays_tot_6m,
  notbaddays_phys_6m = 30 - baddays_phys_6m,
  notbaddays_ment_6m = 30 - baddays_ment_6m)

other_variable_df <- other_variable_df %>% 
  mutate(older = birthyear_list >= 1945 & birthyear_list <= 1958,
         younger = birthyear_list >= 1959 & birthyear_list <= 1989,
         female_cb = female_list==1 & birthyear_list>1968 & birthyear_list<=1989,
         first_week_list = week_list == 1,
         scaled_week_list = week_list - 1)

## 12m timing variables


timing_df <- other_variable_df %>% 
  mutate(dt_returned_0m = round(dt_returned_0m),
         dt_returned_6m = round(dt_returned_6m),
         dt_returned_12m = round(dt_returned_12m)) %>% 
  mutate(mail_to_response_12m = dt_returned_12m - dt_mail_12m) 


# Missing values
timing_df <- timing_df %>% 
  group_by(draw) %>% 
  mutate(mail_to_response_12_mean = mean(ifelse(mail_to_response_12m >= 0, mail_to_response_12m, NA))) %>% 
  mutate(maxmean = round(max(mail_to_response_12_mean))) %>% 
  mutate(mail_to_response_12m = ifelse(mail_to_response_12m <=0, maxmean, mail_to_response_12m)) %>% 
  ungroup() %>% 
  select(-maxmean, -mail_to_response_12_mean) %>% 
  mutate(sample_hdd = 1,
         zip_hh_inc_list = 0)


## Recoding Survey variables


survey_vars_df <- timing_df %>% 
  mutate(health_poor_12m = ifelse(!is.na(health_gen_12m), (health_gen_12m == 1), NA),
         flp_categ_12m = case_when(
           hhinc_pctfpl_12m <= 50 ~ "below 50% FPL",
           between(hhinc_pctfpl_12m, 50, 75) ~ "50-75% FPL",
           between(hhinc_pctfpl_12m, 75, 100) ~ "75-100% FPL",
           between(hhinc_pctfpl_12m, 100, 150) ~ "100-150% FPL",
           150 < hhinc_pctfpl_12m ~ "above 150% FPL"
         ),
         health_genflip_bin_12m = case_when(
           health_gen_bin_12m == 1 ~ 0,
           health_gen_bin_12m == 0 ~ 1
         ),
         health_chgflip_bin_12m = case_when(
           health_chg_bin_12m == 1 ~ 0,
           health_chg_bin_12m == 0 ~ 1
         ),
         health_notpoor_12m = case_when(
           health_poor_12m == 1 ~ 0,
           health_poor_12m == 0 ~ 1
         ),
         notbaddays_tot_12m = 30 - baddays_tot_12m,
         notbaddays_phys_12m = 30 - baddays_phys_12m,
         notbaddays_ment_12m = 30 - baddays_ment_12m)

# happiness
recode_df <- survey_vars_df %>% 
  mutate(happiness_bin_12m = case_when(
    happiness_12m == 1 | happiness_12m == 2 ~ "very/pretty happy",
    happiness_12m == 3 ~ "not too happy"
  ),
  poshappiness_bin_12m = case_when(
    happiness_bin_12m == "very/pretty happy" ~ 0,
    happiness_bin_12m == "not too happy" ~ 1
  ),
  chl_chk_bin_12m = case_when(
    chl_chk_12m == 1 | chl_chk_12m == 2 ~ 1,
    chl_chk_12m == 3 ~ 0
  ),
  dia_chk_bin_12m = case_when(
    dia_chk_12m == 1 | dia_chk_12m == 2 ~ 1,
    dia_chk_12m == 3 ~ 0
  ),
  mam_chk_bin_12m = case_when(
    mam_chk_12m == 1 ~ 1,
    mam_chk_12m == 2 | mam_chk_12m == 3 ~ 0
  ),
  mam_chk_bin_12m = ifelse(birthyear_list > 1968, NA, mam_chk_bin_12m),
  pap_chk_bin_12m = case_when(
    pap_chk_12m == 1 ~ 1,
    pap_chk_12m == 2 | pap_chk_12m == 3 ~ 0
  ))

recode_health_behaviour_df <- recode_df %>% 
  mutate(smk_curr_bin_12m = case_when(
    smk_curr_12m == 1 | smk_curr_12m == 2 ~ 1,
    smk_curr_12m == 3 ~ 0
  ),
  nonsmk_curr_12m = case_when(
    smk_curr_bin_12m == 1 ~ 0,
    smk_curr_bin_12m == 0 ~ 1
  ),
  physical_act_bin_12m = case_when(
    physical_act_12m  == 1 | physical_act_12m == 2 ~ 0,
    physical_act_12m == 3 ~ 1
  ),
  more_active_12m = case_when(
    physical_act_bin_12m == 1 ~ 0,
    physical_act_bin_12m == 0 ~ 1
  ),
  not_er_noner_0m = case_when(
    er_noner_0m == 0 ~ 1,
    er_noner_0m == 1 ~ 0
  ),
  not_er_noner_6m = case_when(
    er_noner_6m == 0 ~ 1,
    er_noner_6m == 1 ~ 0
  ),
  not_er_noner_12m = case_when(
    er_noner_12m == 0 ~ 1,
    er_noner_12m == 1 ~ 0
  ),
  dep_screen_6m = (dep_interest_6m + dep_sad_6m) >= 5,
  dep_screen_12m = (dep_interest_12m + dep_sad_12m) >= 5,
  nodep_screen_12m = case_when(
    dep_screen_12m == 1 ~ 0,
    dep_screen_12m == 0 ~ 1))



cost_tot_oop_99 <- recode_health_behaviour_df %>% 
  select(cost_tot_oop_12m) %>% 
  pull() %>% 
  quantile(0.99, na.rm = TRUE)


cost_tot_owe_99 <- recode_health_behaviour_df %>% 
  select(cost_tot_owe_12m) %>% 
  pull() %>% 
  quantile(0.99, na.rm = TRUE)


recode_health_behaviour_df <- recode_health_behaviour_df %>% 
  mutate(cost_tot_oop_mod_12m = ifelse(cost_tot_oop_12m <= cost_tot_oop_99*2 & cost_tot_oop_12m >= 0, cost_tot_oop_12m, NA),
         cost_tot_owe_mod_12m = ifelse(cost_tot_owe_12m <= cost_tot_owe_99*2 & cost_tot_owe_12m >= 0, cost_tot_owe_12m, NA),
         hhinc_mid_12m = hhinc_cat_12m, 
         hhinc_mid_12m = ifelse(hhinc_mid_12m == 1, 0, hhinc_mid_12m),
         hhinc_mid_12m = ifelse(hhinc_mid_12m != 0, 1250 + 2500 * (hhinc_mid_12m - 2), hhinc_mid_12m),
         hhinc_mid_12m = ifelse(hhinc_cat_12m == 22, 5001, hhinc_mid_12m))

# death variable

final_df <- recode_health_behaviour_df %>% 
  mutate(postn_alive = case_when(
    postn_death == 1 ~ 0,
    postn_death == 0 ~ 1
))

#rm(list = setdiff(ls(), c("final_df", "write_to_disk")))
options(na.action = "na.omit") # Changing back to default NA setting for most people.
if (write_to_disk){
  write.csv(final_df, file = file.path(data_dir, "data_for_analysis_R.csv"), row.names = FALSE)
}


print("Success!")