
results_long <- sim_env$grid_list$results_long
sig_eps_grid <- unique(sim_env$grid_list$results_long$sig_eps)
sig_x_grid <- unique(sim_env$grid_list$results_long$sig_x)

CleanBreaks <- function(grid, n, digits=1) {
    clean_grid <- pretty(grid, n=n)
    clean_grid[1] <- min(grid) + 0.001
    clean_grid <- round(clean_grid, digits=digits)
    return(clean_grid)
}


plot1 <- ggplot(results_long %>% filter(metric == "prop_drop")) +
    geom_tile(aes(x=sig_eps, y=sig_x, fill=value)) +
    geom_tile(aes(x=sig_eps, y=sig_x), fill="gray",
              data=results_long %>% filter(metric == "prop_drop") %>%
                filter(is.na(value))) +
    scale_fill_gradient2(
        trans="log10",
        midpoint=log10(0.01),
        #name=TeX("Estimated $100 \\alpha^*$ to change sign"),
        name=TeX("Percent dropped"),
        breaks=c(0.002, 0.01, 0.05, 0.20),
        labels=c("0.2%", "1%", "5%", "20%")) +
    scale_x_continuous(
      trans="log", breaks=CleanBreaks(sig_eps_grid, n=5), expand=c(0,0)) +
    scale_y_continuous(
      trans="log", breaks=CleanBreaks(sig_x_grid, n=5), expand=c(0,0)) +
    geom_point(aes(x=sig_eps, y=sig_x, shape="bad_value"), size=4,
               data=results_long %>%
                   filter(metric == "prop_drop", value > 0.1)) +
    xlab(TeX("$\\sigma_{\\epsilon}$ (log scale)")) +
    ylab(TeX("$\\sigma_x$  (log scale)")) +
    geom_text(aes(x=min(sig_eps_grid), y=max(sig_x_grid)), label="NA") +
    scale_shape_manual(breaks=c("bad_value"),
                       labels=c("Linearity assumption is suspect"),
                       name=NULL,
                       guide = 'none',
                       values=c(4)) +
    theme(legend.position="bottom", legend.box="vertical") +
    geom_text(aes(x=25, y=4.5),
              label="Linearity assumption suspect\nfor changes > 10%")


accuracy_results_df <- sim_env$acc_list$accuracy_results_df
plot2 <- ggplot(accuracy_results_df) +
   geom_line(
     aes(x=100 * prop_drop, y=diff_pred - diff_true, color="err")) +
   geom_line(aes(x=100 * prop_drop, y=diff_pred, color="lin")) +
   geom_line(aes(x=100 * prop_drop, y=diff_true, color="refit")) +
   xlab("Percent dropped") + ylab("Estimate") +
   theme(legend.position = "bottom") +
   scale_color_discrete(
     name=NULL,
     breaks=c("refit", "lin", "err"),
     labels=c("Actual change", "Linear approximation", "Error"))

grid.arrange(plot1, plot2, ncol=2)
