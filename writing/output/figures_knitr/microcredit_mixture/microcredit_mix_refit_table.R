
refit_df <-
    mcmix_env$results_df %>%
    filter(str_detect(param_name, "^tau"), analysis=="rerun") %>%
    mutate(study_case=droplevels(study_case))

base_df <-
    mcmix_env$base_df %>%
    filter(str_detect(param_name, "^tau")) %>%
    select(analysis, param_name, study_case, metric, value) %>%
    unique() %>%
    mutate(study_case=droplevels(study_case))

table_df <- FormatRefitTable(refit_df, base_df)

RenderLatexTable(
    table_df,
    label="mcmix_re_run_table",
    study_case_label="Model parameter",
    caption=paste0(
        "Microcredit mixture results for a selected set of model parameters. ",
        "Standard errors and ``significance'' are based on the ",
        "estimated 95\\% posterior credible intervals. ",
        # "and stars indicate that zero does not lie in the credible interval.  ",
        GetTableCaptionBoilerplate()
        # "The ``Refit estimate'' column shows the result of re-fitting ",
        # "the model removing ",
        # "the Approximate Most Influential Set. ",
        # "Refits that achieved the desired change are bolded."
    )
)
