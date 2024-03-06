#!/usr/bin/env Rscript

library(tidyverse)
library(latex2exp)
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

paper_directory <- file.path(git_repo_loc, "writing/output/")
data_path <- file.path(paper_directory, "applications_data")
save_path <- file.path(data_path, "simulations")

set.seed(42)

num_obs <- 5000
eps_raw <- rnorm(num_obs)
x_raw <- rnorm(num_obs)
infl_raw <- eps_raw * x_raw / (sum(x_raw^2))
infl_order <- order(infl_raw)
num_obs * infl_raw[infl_order] %>% head()

ProcessInfluence <- function(infl) {
    infl_pos <- infl > 0
    infl_neg <- infl < 0
    
    inds_pos <- (1:length(infl))[infl_pos]
    inds_neg <- (1:length(infl))[infl_neg]
    
    ordered_inds_pos <- inds_pos[order(-1 * infl[infl_pos])]
    ordered_inds_neg <- inds_neg[order(infl[infl_neg])]
    
    infl_cumsum_pos <- cumsum(infl[ordered_inds_pos])
    infl_cumsum_neg <- cumsum(infl[ordered_inds_neg])
    return(list(
        neg=list(infl_inds=ordered_inds_neg,
                 infl_cumsum=infl_cumsum_neg),
        pos=list(infl_inds=ordered_inds_pos,
                 infl_cumsum=infl_cumsum_pos),
        num_obs=length(infl)
    ))
}



GetAlpha <- function(infl_dat, delta) {
    # To produce a negative change, drop observations with positive influence
    # scores, and vice-versa.
    infl_list <- if (delta < 0) infl_dat$pos else infl_dat$neg
    num_obs <- length(infl_list$infl_cumsum)
    n_vec <- 1:num_obs
    n_drop <- approx(x=-1 * c(0, infl_list$infl_cumsum), 
                     y=c(0, n_vec), 
                     xout=delta)$y %>% ceiling()
    if (is.na(n_drop)) {
        drop_inds <- NA
    } else if (n_drop == 0) {
        drop_inds <- c()
    } else {
        drop_inds <- infl_list$infl_inds[1:n_drop]
    }
    return(list(
        n=n_drop,
        prop=n_drop / infl_dat$num_obs,
        inds=drop_inds
    ))
}


if (FALSE) {
  # Test
  infl_test <- c(-3, 5, -2, 8, 0, -1)
  infl_dat <- ProcessInfluence(infl_test)
  infl_dat
  
  infl_test[infl_dat$neg$infl_inds]
  infl_test[infl_dat$pos$infl_inds]
  
  GetAlpha(infl_dat, 1)
  GetAlpha(infl_dat, 4)
  GetAlpha(infl_dat, 5.5)
  GetAlpha(infl_dat, 7)
  
  GetAlpha(infl_dat, -1)
  GetAlpha(infl_dat, -12)
  GetAlpha(infl_dat, 14)
  GetAlpha(infl_dat, 0)
}


ScaleData <- function(infl_raw, x_raw, eps_raw, theta0, sig_x, sig_eps) {
    x <- sig_x * x_raw
    eps <- sig_eps * eps_raw
    y <- x * theta0 + eps
    infl <- eps * x / (sum(x^2))
    return(list(
        num_obs=length(eps_raw),
        infl=infl,
        infl_dat=ProcessInfluence(infl),
        x=x,
        y=y,
        theta0=theta0,
        thetahat=sum(y * x) / sum(x * x)))
}

DropRows <- function(data_list, drop_inds) {
    diff_pred <- -1 * sum(data_list$infl[drop_inds])
    keep_inds <- setdiff(1:data_list$num_obs, drop_inds)
    diff_true <-
        sum(data_list$x[keep_inds] * data_list$y[keep_inds]) / sum(data_list$x[keep_inds]^2) -
        data_list$thetahat
        
    return(list(diff_pred=diff_pred,
                diff_true=diff_true,
                thetahat=data_list$thetahat,
                theta_pred=diff_pred + data_list$thetahat,
                theta_true=diff_true + data_list$thetahat))
}

#######################################

grid_points <- 20
sig_x_grid <- exp(seq(0.1, 2, length.out=grid_points))
sig_eps_grid <- exp(seq(2, 4, length.out=grid_points))

theta0 <- 0.5
alpha <- 0.01

ScaleThisData <- function(theta0, sig_x, sig_eps) {
    return(ScaleData(infl_raw=infl_raw, x_raw=x_raw, eps_raw=eps_raw,
                     theta0=theta0, sig_x=sig_x, sig_eps=sig_eps))
}

results_df <- data.frame()
pb <- txtProgressBar(min=0, max=length(sig_x_grid) * length(sig_eps_grid), style=3)
i <- 1
for (sig_x in sig_x_grid) { for (sig_eps in sig_eps_grid) {
    setTxtProgressBar(pb, i)
    data_list <- ScaleThisData(theta0, sig_x=sig_x, sig_eps=sig_eps)
    amis_list <- GetAlpha(data_list$infl_dat, data_list$thetahat)
    if (!is.na(amis_list$n)) {
        drop_result <- DropRows(data_list, amis_list$inds) %>% data.frame()
    } else {
        drop_result <- data.frame(diff_pred=NA)
    }
    i <- i + 1
    results_df <- bind_rows(
        results_df,
        drop_result %>%
            mutate(theta0=theta0, 
                   sig_x=sig_x, 
                   sig_eps=sig_eps,
                   noise=sig_eps^2 / sig_x^2,
                   thetahat=data_list$thetahat,
                   num_drop=amis_list$n,
                   prop_drop=amis_list$prop)
        )
}}
close(pb)


#summary(results_df$thetahat)
#names(results_df)

results_long <- pivot_longer(
    results_df, 
    cols=-c(sig_x, sig_eps, noise, thetahat), names_to="metric")

#results_long %>% filter(metric == "prop_drop") %>% pull("value") %>% max()

grid_list <- list()
grid_list[["results_long"]] <- results_long
grid_list[["results_df"]] <- results_df
grid_list[["sig_x_grid"]] <- sig_x_grid
grid_list[["sig_eps_grid"]] <- sig_eps_grid
grid_list[["num_obs"]] <- num_obs
grid_list[["theta0"]] <- theta0


if (FALSE) {
    CleanBreaks <- function(grid, n, digits=1) {
        clean_grid <- pretty(grid, n=n)
        clean_grid[1] <- min(grid) + 0.001
        clean_grid <- round(clean_grid, digits=digits)
        return(clean_grid)
    }

    ggplot(results_long %>% filter(metric == "prop_drop")) +
    geom_tile(aes(x=sig_eps, y=sig_x, fill=value)) +
    geom_tile(aes(x=sig_eps, y=sig_x), fill="gray",
              data=results_long %>% filter(metric == "prop_drop") %>% filter(is.na(value))) +
    scale_fill_gradient2(
        trans="log10",
        midpoint=log10(0.01),
        name=TeX("Estimated $\\alpha^*$ to change sign")) +
    scale_x_continuous(trans="log", breaks=CleanBreaks(sig_eps_grid, n=5), expand=c(0,0)) +
    scale_y_continuous(trans="log", breaks=CleanBreaks(sig_x_grid, n=5), expand=c(0,0)) +
    geom_point(aes(x=sig_eps, y=sig_x, shape="bad_value"), size=4,
               data=results_long %>%
                   filter(metric == "prop_drop", value > 0.1)) +
    xlab(TeX("$\\sigma_{\\epsilon}$ (log scale)")) +
    ylab(TeX("$\\sigma_x$  (log scale)")) +    
    geom_text(aes(x=min(sig_eps_grid), y=max(sig_x_grid)), label="NA") +
    scale_shape_manual(breaks=c("bad_value"),
                       labels=c("Linearity assumption is suspect"),
                       name=NULL,
                       values=c(4))

}


#######################################

sig_x <- 2
sig_eps <- 1

theta0 <- 0
data_list <- ScaleThisData(0, sig_x=sig_x, sig_eps=sig_eps)
alpha_max <- 0.1
alpha_grid <- exp(seq(log(1 / data_list$num_obs), log(alpha_max), length.out=500))
n_grid <- ceiling(data_list$num_obs * alpha_grid) %>% unique()

accuracy_results_df <- data.frame()
pb <- txtProgressBar(min=0, max=length(n_grid), style=3)
i <- 1
for (n_drop in n_grid) {
    setTxtProgressBar(pb, i)
    drop_result <- DropRows(
        data_list, 
        data_list$infl_dat$neg$infl_inds[1:n_drop]) %>% 
        data.frame()
    i <- i + 1
    accuracy_results_df <- bind_rows(
        accuracy_results_df,
        drop_result %>%
            mutate(num_drop=n_drop,
                   prop_drop=n_drop / data_list$num_obs)
    )
}
close(pb)


if (FALSE) {
    ggplot(accuracy_results_df) +
        geom_line(aes(x=prop_drop, y=diff_pred - diff_true, color="Linear Approximation Error")) +
        geom_line(aes(x=prop_drop, y=diff_pred, color="Linear Approximation")) +
        geom_line(aes(x=prop_drop, y=diff_true, color="Actual change")) +
        xlab(TeX("$\\alpha$")) + ylab("")
}

acc_list <- list()
acc_list[["theta0"]] <- theta0
acc_list[["num_obs"]] <- num_obs
acc_list[["sig_x"]] <- sig_x
acc_list[["alpha_max"]] <- alpha_max
acc_list[["sig_eps"]] <- sig_eps
acc_list[["accuracy_results_df"]] <- accuracy_results_df



#######################################
# Save

save(grid_list, acc_list,
     file=file.path(save_path, "noise_grid.Rdata"))

print("Success!")