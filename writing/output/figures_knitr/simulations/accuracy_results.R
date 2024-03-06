
accuracy_results_df <- sim_env$acc_list$accuracy_results_df
ggplot(accuracy_results_df) +
    geom_line(aes(x=prop_drop, y=diff_pred - diff_true,
              color="Linear Approximation Error")) +
    geom_line(aes(x=prop_drop, y=diff_pred, color="Linear Approximation")) +
    geom_line(aes(x=prop_drop, y=diff_true, color="Actual change")) +
    xlab(TeX("$\\alpha$")) + ylab("")
