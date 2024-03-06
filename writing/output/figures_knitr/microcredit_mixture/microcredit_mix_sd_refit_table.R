# Ugh, make a custom table for the log sd tau parameters

fontsize <- "\\scriptsize"
label <- "mcmix_sd_re_run_table"
study_case_label <- "Model parameter"

caption <- paste0(
    "Results for the log posterior standard deviation estimates of ",
    "the effect size distribution in the microcredit mixture model.  ",
    "Sign and significance are not meaningful for ",
    "posterior standard deviations, so we drop 0.5\\% of datapoints to ",
    "attempt to produce large positive and negative changes.  ",
    "The ``Refit estimate'' column shows the result of re-fitting ",
    "the model removing ",
    "Approximate Most Influential Set.  ",
    "The ``Prediction'' column shows the predicted change under ",
    "the same perturbation.  "
)

results_df <-
    mcmix_env$results_df %>%
    filter(str_detect(param_name, "^log_sd_tau")) %>%
    filter(description %in% sprintf("%s177", c("neg", "pos"))) %>%
    mutate(study_case=droplevels(study_case))

base_df <-
    mcmix_env$base_df %>%
    filter(str_detect(param_name, "^log_sd_tau")) %>%
    select(analysis, param_name, study_case, metric, value) %>%
    unique() %>%
    mutate(study_case=droplevels(study_case))


# Rename the change_string for the custom change descriptions.
# Note that "neg" and "pos" refer to the sign of the influence
# function, which the opposite of the direction produced when
# points are dropped.
table_df <-
    FormatRefitTable(results_df %>% filter(analysis %in% c("rerun", "pred")),
                     base_df %>% select(-analysis)) %>%
    mutate(change_string=ordered(
        description,
        levels=c("neg177", "pos177"),
        labels=c("Increase",
                 "Decrease")))
    # mutate(change_string=ordered(
    #     description,
    #     levels=c("neg177", "pos177"),
    #     labels=c("Drop 0.5\\% to increase",
    #              "Drop 0.5\\% to decrease")))



# A custom version of RenderNewLatexTable
cat("\n\n")
cat("\\begin{table}")
cat(fontsize, "\n")
cat("\\begin{tabular}{|cccccc|}\n")
cat("\\toprule\n")
cat(study_case_label,
    " & Original ",
    " & Change ",
    " & Estimate  ",
    " & Prediction ",
    " & Observations ",
    "\\\\\n", sep="")
cat(" ",
    " & estimate ",
    " & type ",
    " & after refit ",
    " &  ",
    " & dropped ",
    "\\\\\n", sep="")
cat("\\midrule\n")
cat("\\midrule\n")


# cat("\n\n")
# cat("\\begin{table}")
# cat(fontsize, "\n")
# cat("\\begin{tabular}{|cccccc|}\n")
# cat("\\toprule\n")
# cat(study_case_label,
#     " & Original estimate ",
#     " & Change type ",
#     " & Refit estimate ",
#     " & Prediction ",
#     " & Observations dropped ",
#     "\\\\\n", sep="")
# cat("\\midrule\n")
# cat("\\midrule\n")


for (case in levels(table_df$study_case)) {
    this_table_df <- filter(table_df, study_case == case)
    refit_df <- filter(this_table_df, analysis == "rerun")
    pred_df <- filter(this_table_df, analysis == "pred")

    # All the estimate strings for the original should be the same
    stopifnot(length(unique(this_table_df$est_string_base)) == 1)
    est_string_base <- sprintf("%0.3f", this_table_df$param_base[1])

    # There should be two target changes in the refit
    stopifnot(nrow(refit_df) == 2)
    stopifnot(nrow(pred_df) == 2)

    for (i in 1:2) {
        if (i == 1) {
            cat(case, " & ", est_string_base)
        } else {
            cat("&")
        }
        stopifnot(refit_df[["change_string"]][i] == pred_df[["change_string"]][i])
        param_refit <- sprintf("%0.3f", refit_df[["param_rerun"]][i])
        param_pred <- sprintf("%0.3f", pred_df[["param_rerun"]][i])

        cat(" & ", refit_df[["change_string"]][i] %>% as.character())
        # cat(" & ", refit_df[["param_rerun"]][i] %>% as.character())
        # cat(" & ", pred_df[["param_rerun"]][i] %>% as.character())
        cat(" & ", param_refit)
        cat(" & ", param_pred)
        cat(" & ", refit_df[["num_string"]][i] %>% as.character())
        cat("\\\\\n")
    }
    cat("\\midrule\n")
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
if (!is.null(caption)) {
    cat("\\caption{", caption, "}\n")
}
if (!is.null(label)) {
    cat("\\label{table:", label, "}\n", sep="")
}
cat("\\end{table}")

cat("\n\n")
