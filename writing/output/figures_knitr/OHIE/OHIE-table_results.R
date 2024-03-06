
outcome_df <- bind_rows(
    data.frame(outcome="health_genflip_bin_12m", outcome_string="Health genflip"),
    data.frame(outcome="health_notpoor_12m", outcome_string="Health notpoor"),
    data.frame(outcome="health_chgflip_bin_12m", outcome_string="Health change flip"),
    data.frame(outcome="notbaddays_tot_12m", outcome_string="Not bad days total"),
    data.frame(outcome="notbaddays_phys_12m", outcome_string="Not bad days physical"),
    data.frame(outcome="notbaddays_ment_12m", outcome_string="Not bad days mental"),
    data.frame(outcome="nodep_screen_12m", outcome_string="Nodep Screen"),
)
stopifnot(setequal(outcome_df$outcome, unique(ohie_env$results_df$outcome)))

AddStudyCase <- function(df) {
  df %>%
    mutate(
        outcome_string=ordered(
            outcome,
            levels=outcome_df$outcome,
            labels=outcome_df$outcome_string),
        study_case=outcome_string) %>%
    return()
}

table_df <- AddStudyCase(
  ohie_env$results_df %>% filter(analysis=="rerun"))
base_df <- AddStudyCase(ohie_env$base_df)

filter(table_df, method == "iv") %>%
  FormatRefitTable(base_df) %>%
  RenderLatexTable(
    label="ohie_profit_results_iv",
    caption=paste0(
        "Medicaid profit results with IV for a range of outcome ",
        "variables.  ",
        GetTableCaptionBoilerplate())
    )


filter(table_df, method == "regression") %>%
  FormatRefitTable(base_df) %>%
  RenderLatexTable(
    label="ohie_profit_results_reg",
    caption=paste0(
        "Medicaid profit results with OLS for a range of outcome ",
        "variables.  ",
        "All outcomes are twelve-month outcomes.  ",
        GetTableCaptionBoilerplate())
    )
