library(tidyverse)
library(gridExtra)
library(latex2exp)
library(zaminfluence)

base_dir <- "/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench"
#base_dir <- "/Users/rachaelmeager/AdversarialInfluenceWorkbench"

paper_directory <- file.path(base_dir, "writing/output/")
data_path <- file.path(paper_directory, "applications_data")

# Save an Rdata file in save_path
save_path <- file.path(data_path, "simulations")
save_list <- list()

# Optionally save static figures in output_dir
save_static_figures <- FALSE
output_dir <- file.path(base_dir, "writing/output/static_figures")



GetDf <- function(x, raw_xi, sigma) {
    gamma <- raw_xi / sqrt(mean(raw_xi^2))
    xi <- sigma * gamma
    df <- data.frame(xi=xi, gamma=gamma) %>%
        pivot_longer(cols=c(xi, gamma)) %>%
        group_by(name) %>%
        mutate(i=order(value))
    return(df)
}

n <- 20
# x <- rnorm(floor(n / 2))
# x <- c(-x, x) %>% sort()
x <- 1:20
x <- x - mean(x)
x <- x / sd(x)

ymin <- -5
ymax <- 5

# l, h = light, heavy
# w, n = wide, narrow

GetHeavyXi <- function(x) {
    #return(sign(x) * x ^ 2)
    return(sign(x) / (abs(x) ^ (1.5)))
}

GetLightXi <- function(x) {
    return(pnorm(x) - 0.5)
}

# Set some sigma to remain within the limits of the graph.
light_xi <- GetLightXi(x)
heavy_xi <- GetHeavyXi(x)
sig_w <- 4.5 / max(heavy_xi / sqrt(mean(heavy_xi^2)))
sig_n <- 0.5 / max(light_xi / sqrt(mean(light_xi^2)))

PlotInfluence1 <- function(df) {
    ggplot(df) +
        geom_bar(aes(x=i, y=value, fill=name),
                 stat="identity", position="dodge", color="black") +
        scale_fill_grey() +
        ylim(ymin, ymax)
}

PlotInfluence1(GetDf(x, GetHeavyXi(x), sig_w))

df <- GetDf(x, GetHeavyXi(x), sig_n)
g_hn <- PlotInfluence1(df)

df <- GetDf(x, GetHeavyXi(x), sig_w)
g_hw <- PlotInfluence1(df)


df <- GetDf(x, GetLightXi(x), sig_n)
g_ln <- PlotInfluence1(df)

df <- GetDf(x, GetLightXi(x), sig_w)
g_lw <- PlotInfluence1(df)

grid.arrange(
    g_hn, g_hw, g_ln, g_lw, ncol=2
)




##########################

PlotInfluence2 <- function(df, i_max) {
    df_plot <- df %>% filter(name == "xi")
    ggplot(df_plot) +
        geom_bar(aes(x=i, y=value), stat="identity") +
        geom_bar(aes(x=i, y=value), stat="identity",
                 data=filter(df_plot, i <= i_max), fill="red") +
        ylim(ymin, ymax) +
        xlab("Influence rank") + ylab("Influence")
}

i_max <-4

df_hw <- GetDf(x, GetHeavyXi(x), sig_w)
g_hw_effect <- -sum(filter(df_hw, name == "xi", i <= i_max)$value)  / n
g_hw <-
    PlotInfluence2(df_hw, i_max) +
    ggtitle(TeX(sprintf("Heavy tails. $\\Gamma_{\\alpha}$ = %0.3f", g_hw_effect))) +
    ylab(TeX("$\\gamma_n$"))

df_lw<- GetDf(x, GetLightXi(x), sig_w)
g_lw_effect <- -sum(filter(df_lw, name == "xi", i <= i_max)$value) / n
g_lw <-
    PlotInfluence2(df_lw, i_max) +
    ggtitle(TeX(sprintf("Light tails. $\\Gamma_{\\alpha}$ = %0.3f", g_lw_effect))) +
    ylab(TeX("$\\gamma_n$"))

if (save_static_figures) {
    png(width=6, height=4, res=300, units="in", filename=file.path(output_dir, "gamma_alpha_example.png"))
    grid.arrange(g_hw, g_lw, ncol=2)
    dev.off()
}

##########################


PlotInfluence3 <- function(df, x_max, xlim=4, bins=50, delta=-0.2, scale=1.0, title=NULL) {
    df_plot <- df %>% mutate(x=x * scale)
    drop_rows <- df_plot$x <= x_max * scale
    infl <- -1 * sum(df_plot[drop_rows, "x"] + delta) / nrow(df_plot)
    alpha <- mean(drop_rows)
    if (is.null(title)) {
        title <- sprintf("N = %d   Influence = %0.3f   Alpha = %0.3f", nrow(df), infl, alpha)
    }
    ggplot(df_plot) +
        geom_histogram(aes(x=x + delta, y=..count..), bins=bins) +
        geom_histogram(aes(x=x + delta, y=..count..), bins=bins,
                 data=filter(df_plot, x <= x_max * scale), fill="red") +
        xlim(-xlim, xlim) +
        geom_vline(aes(xintercept=delta + infl), color="red") +
        geom_vline(aes(xintercept=delta), color="blue") +
        geom_vline(aes(xintercept=0), color="black") +
        xlab(TeX("$\\epsilon$")) + ylab("") +
        ggtitle(title)
}

set.seed(50)
dfs <- data.frame(x=rnorm(100))
dfm <- data.frame(x=rnorm(1000))
dfl <- data.frame(x=rnorm(10000))

dfs_sort <- sort(dfs$x)
dfm_sort <- sort(dfm$x)
dfl_sort <- sort(dfl$x)

if (save_static_figures) {
    PlotInfluence3(dfs, -1e6, bins=30, title="")
    ggsave("~/Downloads/base_plot.png", width=5, height=4, units="in")
    
    PlotInfluence3(dfs, dfs_sort[1], bins=30, title="")
    ggsave("~/Downloads/base_plot_remove1.png", width=5, height=4, units="in")
    
    PlotInfluence3(dfs, dfs_sort[10], bins=30, title="")
    ggsave("~/Downloads/base_plot_remove10.png", width=5, height=4, units="in")
    
    PlotInfluence3(dfl, dfl_sort[1000], bins=200)
    ggsave("~/Downloads/dfl_plot_remove0.1.png", width=5, height=4, units="in")
}

if (save_static_figures) {
    png(width=5, height=4, res=300, units="in",
        filename=file.path("~/Downloads/varying_n_plot_remove1.png"))
    grid.arrange(
        PlotInfluence3(dfs, dfs_sort[1], bins=100),
        PlotInfluence3(dfm, dfm_sort[1], bins=150),
        PlotInfluence3(dfl, dfl_sort[1], bins=200),
        ncol=1)
    dev.off()
    
    png(width=5, height=4, res=300, units="in", filename=file.path("~/Downloads/varying_n_plot_remove10.png"))
    grid.arrange(
        PlotInfluence3(dfs, dfs_sort[10], bins=100),
        PlotInfluence3(dfm, dfm_sort[10], bins=150),
        PlotInfluence3(dfl, dfl_sort[10], bins=200),
        ncol=1)
    dev.off()
    
    
    png(width=5, height=4, res=300, units="in", filename=file.path("~/Downloads/varying_n_plot_remove0.01.png"))
    grid.arrange(
        PlotInfluence3(dfs, dfs_sort[10], bins=200),
        PlotInfluence3(dfm, dfm_sort[100], bins=200),
        PlotInfluence3(dfl, dfl_sort[1000], bins=200),
        ncol=1)
    dev.off()
}



grid.arrange(
grid.arrange(
    PlotInfluence3(dfs, dfs_sort[1]),
    PlotInfluence3(dfm, dfm_sort[1]),
    PlotInfluence3(dfl, dfl_sort[1]),
    ncol=1)
,
grid.arrange(
    PlotInfluence3(dfs, dfs_sort[10]),
    PlotInfluence3(dfm, dfm_sort[10]),
    PlotInfluence3(dfl, dfl_sort[10]),
    ncol=1)
,
grid.arrange(
    PlotInfluence3(dfs, dfs_sort[10]),
    PlotInfluence3(dfm, dfm_sort[100]),
    PlotInfluence3(dfl, dfl_sort[1000]),
    ncol=1)
, ncol=3
)

if (save_static_figures) {
    png(width=5, height=4, res=300, units="in", filename=file.path("~/Downloads/varying_scale_plot.png"))
    grid.arrange(
        PlotInfluence3(dfl, dfl_sort[1000], scale=0.5, bins=400),
        PlotInfluence3(dfl, dfl_sort[1000], bins=400),
        PlotInfluence3(dfl, dfl_sort[1000], scale=1.5, bins=400),
        ncol=1)
    dev.off()
}



dfunif <- data.frame(x=runif(10000) - 0.5)
Renorm <- function(x1, x2) {
    x1 <- x1 / sqrt(sum(x1^2))
    return(x1 * sqrt(sum(x2^2)))
}

dfunif2 <- dfunif
dfunif2$x <- sign(dfunif2$x) * sqrt(abs(dfunif2$x))

dfunif$x <- Renorm(dfunif$x, dfl$x)
dfunif2$x <- Renorm(dfunif2$x, dfl$x)

dfunif_sort <- sort(dfunif$x)
dfunif2_sort <- sort(dfunif2$x)

dfexp <- data.frame(x=-1 * rexp(10000))
dfexp$x <- dfexp$x - mean(dfexp$x)
dfexp$x <- Renorm(dfexp$x, dfl$x)
dfexp_sort <- sort(dfexp$x)

if (save_static_figures) {
    png(width=5, height=4, res=300, units="in", filename=file.path("~/Downloads/varying_shape_plot.png"))
    grid.arrange(
        PlotInfluence3(dfl, dfl_sort[1000], bins=400),
        PlotInfluence3(dfunif2, dfunif2_sort[1000], bins=400),
        PlotInfluence3(dfexp, dfexp_sort[1000], bins=400),
        ncol=1)
    dev.off()
}



###############################
###############################
###############################
# Super simple plot.  This is the one we actually use.

TransformXi <- function(x) {
    abs_x <- abs(x)
    return(sign(x) * abs_x^(1.15))
}

set.seed(50)
df_plot <- data.frame(x=sort(TransformXi(rnorm(10000))))
alpha <- 0.05
x_max <- quantile(df_plot$x, alpha)
df_plot <- df_plot %>%
    mutate(drop= x <= x_max,
           alpha=!!alpha)

bins <- 300
xlim <- max(abs(df_plot$x))

ggplot(df_plot) +
    geom_histogram(aes(x=x, y=..count..), bins=bins) +
    geom_histogram(aes(x=x, y=..count..), bins=bins,
                   data=filter(df_plot, drop), fill="red") +
    xlim(-xlim, xlim) +
    xlab(TeX("$\\psi$")) + ylab("") +
    ggtitle(TeX("Influence score histogram (N = 10000, $\\alpha$ = 0.05)"))

save_list[["df"]] <- df_plot

if (save_static_figures) {
    png(width=5, height=3.5, res=300, units="in", filename=file.path("~/Downloads/simple_infl_example.png"))
    ggplot(df_plot) +
        geom_histogram(aes(x=x, y=..count..), bins=bins) +
        geom_histogram(aes(x=x, y=..count..), bins=bins,
                       data=filter(df_plot, x <= x_max), fill="red") +
        xlim(-xlim, xlim) +
        xlab(TeX("$\\psi$")) + ylab("") +
        ggtitle(TeX("Influence score histogram (N = 10000, $\\alpha$ = 0.05)"))
    dev.off()
}



#######################################
# Optionally save an Rdata file

if (FALSE) {
    save(save_list, file=file.path(save_path, "visualization.Rdata"))
}


