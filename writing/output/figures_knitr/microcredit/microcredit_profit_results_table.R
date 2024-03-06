results_df <-
    microcredit_env$results_df  %>%
    filter(analysis=="rerun") %>%
    mutate(study_case=ordered(site))

base_df <-
    microcredit_env$base_df %>%
    mutate(study_case=ordered(site))

FormatRefitTable(results_df, base_df) %>%
  RenderLatexTable(
    label="mc_profit_results",
    caption=paste0(
      "Microcredit regressions for the profit outcome. ",
      GetTableCaptionBoilerplate())
  )
