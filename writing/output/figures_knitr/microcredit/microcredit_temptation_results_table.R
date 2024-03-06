
results_df <-
    microcredit_temptation_env$results_df %>%
    filter(analysis=="rerun") %>%
    mutate(study_case=ordered(site))

base_df <-
    microcredit_temptation_env$base_df %>%
    mutate(study_case=ordered(site))

FormatRefitTable(results_df, base_df) %>%
  RenderLatexTable(
    label="mc_temptation_results",
    caption=paste0(
      "Microcredit regressions for the temptation outcome. ",
      GetTableCaptionBoilerplate())
  )

#
# microcredit_temptation_env$results_df %>%
#     mutate(study_case=ordered(site)) %>%
#     FormatNewRefitTable() %>%
#     RenderNewLatexTable(
#       label="mc_temptation_results",
#       caption=paste0(
#         "Microcredit regressions for the temptation outcome. ",
#         GetTableCaptionBoilerplate())
#     )
