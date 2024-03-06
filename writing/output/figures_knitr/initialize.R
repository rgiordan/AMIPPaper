# Initialize R for knitr.

library(tidyverse)
library(knitr)
library(xtable)
library(gridExtra)
library(latex2exp)
library(zaminfluence)
library(stargazer)

# This must be run from within the git repo, obviously.
git_repo_loc <- system("git rev-parse --show-toplevel", intern=TRUE)
paper_directory <- file.path(git_repo_loc, "writing/output/")
data_path <- file.path(paper_directory, "applications_data")
source(file.path(paper_directory, "figures_knitr", "table_formatting_lib.R"))

# Set some figure defaults.
# opts_chunk$set(fig.width=4.9, fig.height=3)
opts_chunk$set(fig.pos='!h', fig.align='center', dev='png', dpi=300)
opts_chunk$set(echo=knitr_debug, message=knitr_debug, warning=knitr_debug)

# Set the default ggplot theme.
theme_set(theme_bw())

# Load into an environment rather than the global space.
LoadIntoEnvironment <- function(filename) {
  my_env <- environment()
  load(filename, envir=my_env)
  return(my_env)
}

# Width in rendered units.  High numbers give smaller fonts on the images.
# Setting base_image_width to the actual width of a page of paper (8.5)
# lets ggplot pick what it thinks is best.
base_aspect_ratio <- 10 / 8.5 # Basically fill a page by default
base_image_width <- 8.5
#base_image_width <- 7

SetImageSize <- function(aspect_ratio, image_width=base_image_width) {
  # Set the size on the page
  ow <- "0.98\\linewidth"
  oh <- sprintf("%0.3f\\linewidth", aspect_ratio * 0.98)

  # Set the size in rendering
  fw <- image_width
  fh <- image_width * aspect_ratio

  opts_chunk$set(out.width=ow,
                 out.height=oh,
                 fig.width=fw,
                 fig.height=fh)
}


SetFullImageSize <- function() SetImageSize(
    aspect_ratio=base_aspect_ratio, image_width=base_image_width)

# Default to a full image.
SetFullImageSize()

# A convenient funciton for extracting only the legend from a ggplot.
# Taken from
# https://tinyurl.com/y8c742p6
GetLegend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Define common colors.
GGColorHue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


GetGraphColors <- function(legend_breaks) {
  stopifnot(length(legend_breaks) <= 4)
  graph_colors <- GGColorHue(4)[1:length(legend_breaks)]
  names(graph_colors) <- legend_breaks
  return(graph_colors)
}
