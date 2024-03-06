# Knitr instructions

The old `figures.tex` file has been replaced with the knitr file
`figures_knitr.Rnw`.  To use it, run `./knit_to_tex.sh`.  This will
generate `figures_knitr.tex`, and ordinary `tex` file.  With this file in
hand, you can run `pdflatex robustness-draft-RM.tex` as normal.

Do not edit `figures_knitr.tex` directly.  Your changes will be overwritten the
next time you run `./knit_to_tex.sh`.  To make changes to the figures, edit
`figures_knitr.Rnw` or its R include files, and then re-run `./knit_to_tex.sh`.
Maybe the most common mode of failure with knitr is forgetting this, making a
bunch of changes to `figures_knitr.tex`, and then losing them!

## How it works

The knitr file `figures_knitr.Rnw` mostly consists of "code chunks", which begin
with an `<<...>>=` and end wtih an `@`.  Everything oustide these code chunks is
ordinary LaTeX and is unchanged when running knitr.  Running knitr executes the
R code inside the chunks and puts the result into the tex file which it outputs.

I don't like to include R directly in the LaTeX because it's hard to debug,
so you'll see that my code chunks mostly just source `.R` files from the
directory `figures_knitr`.  So to change your graphs, you can mostly
just change the pure R code in the sourced scripts, re-knit, and you will
see the changes in the final results.  (But make sure you have turned
off caching---see below.)

### Plot and table chunks

The knitr chunk settings for plots and tables are different.  I have set
them up correctly, and you can just copy and paste them, but make sure
you get the settings for the correct chunk type.

Tables work using R functions that output TeX code directly.  The
`results='asis'` chunk option simply prints the output of the R code
to the TeX file.

Plots are saved as `png` files in the `figures/` directory, and then imported
using `includegraphics`.  Knitr takes care of all the float environment stuff. I
always use two chunks per plot.  The first sets the caption for the plot,
storing it in the R vairable `figcap`.  This is nice because it avoids having
lengthy strings in the argument to a chunk.  It also allows you to set the
figure size in the caption chunk, as I describe in the next point.

### Setup script: `initialize.R` and figure sizes.

The script `initialize.R` loads all the libraries you'll need to generate
your graphs and tables.  It also defines some helper functions.
The most important on sets figure sizes, and it's called
`SetImageSize`.

LaTeX + R control images with four parameters: rendering height and width, and
printing height and width.  The first controls how big ggplot thinks the plot
window is, and the second controls how the plot is displayed on the page.  You
want the aspect ratio (height / width) of the two to be the same, otherwise the
image will appear distorted.

The function `SetImageSize` allows you to specify two parameters to control
the size of a plot: the ggplot width, and the aspect ratio.  It then
makes sure the aspect ratios match, and that the plot takes up the whole
width of the printed page.

The default aspect ratio and image width are set in `initialize.R`:

```
base_aspect_ratio <- 8 / (5 * 2)

base_image_width <- 8.5
```

You can change these defaults for everything simply by changing these values.
If `base_aspect_ratio` is larger, plots will be taller.  If `base_image_width`
is larger, the fonts will be smaller (because `ggplot` thinks it has a
bigger window to print in).

To change the size for a single plot, insert the command `SetImageSize` into a
chunk _before_ the plot chunk.  The figures size is set when the plot chunk is
called, so it does no good to set the figure size inside the same chunk as where
you are producing your plot. Since I always define a caption chunk for a figure,
you can put `SetImageSize` in the caption chunk.  Note that a call to
`SetImageSize` will set the size for all subsequent figures as well, so you need
to set it back to the default after you change it.  I define the convenient
`SetFullImageSize()` to set the size back to the default.  I'm pretty sure you
can call it within a plotting chunk, after you make the plot, and it will apply
to the next chunk.

The image size controls the size of the entire figure.  If you use
`grid.arrange` (as I do often) containing, say, six sub-figures, the image
size set by knitr will be the size for the entire six-figure plot.  For this
reason, at the default size, a six-figure `grid.arrange` plot and a single
stand-alone plot will appear to be at very different resolutions.  If you
want to avoid this, you can set a custom larger size for plots which you
know will contain multiple plots.  You might even define a `SetImageGridSize`
function that automatically sets an appropriate aspect ratio for a figure
that you know contains `x` rows and `y` columns.

### Stand-alone legends

Sometimes you don't want to repeat a legend in a figure containing
many plots.  You can use the `GetLegend` function, which extracts
from a ggplot object a stand-alone plot containing only the legend.
To use it, generate one plot with a legend, extract it, and then re-generate
that plot and all the others without a legend.  Then place the extracted
legend at the appropriate place in a call to, say, `grid.arrange`.

### Data loading script: `load_data.R`

Each analysis has its data loaded into a separate environment variable.
This is done in `load_data.R`.  When copying over graphs from the original
analyses, make sure you copy the variables you need from the analysis
environment first, or knitr won't be able to find them.

### Debugging and designing: `debugging.R`

The script `debugging.R` does not get called by knitr.  Instead, it's an
R script that sets the variables you need to run knitr, so you can source
your graphs in R and see the output without rendering LaTeX.  I like to
tweak my graphs using `debugging.R` since it's faster and pure R.  The
results automatically show up in LaTeX since knitr sources the same files
you're editing.

### Debugging knitr

If something is going wrong, you can often tell what it is by looking at
the output, `figures_knitr.tex`.  When there are errors, the TeX may not
compile, so maybe start by looking at the source.

By default, knitr doesn't produce much output (because you don't want all
your R commands showing up in your document).  For debugging, however,
you can output all the R stuff to TeX by setting `knitr_debug <- TRUE`
in the `setup` chunk of `figures_knitr.Rnw`.  Don't forget to set it back
to `FALSE` when you're done debugging.

### Caching

To speed up the knit, I set caching variables in the `setup` block:

```
cash_cache <- TRUE
ohie_cache <- TRUE
microcredit_cache <- TRUE
mcmix_cache <- TRUE
sim_cache <- TRUE
```

If an analysis's cache variable is set to true, its graphs and tables are not
regenerated when you knit!  This saves you time when you're only fixing one of
the analyses, but you have to remember to set a cache variable to `FALSE` is you
want to see the effect of changes.  (And set it back to `TRUE` when you're done
if you don't want to re-generate the figures every time.)

You can also always manually clear the cache by running `rm figure/*`.

### R variables

Note that knitr is bascially running an R shell as it processes the knitr
file, and you can use R variable in chunk definitions.  If you find yourself
doing something repetitively (like setting a bunch of captions to have
the same text), note that you could define a single R variable and
re-use it.

### References

You can refer to a chunk by its name and the appropriate prefix.  For
example, you can refer to `cash_transfers_sensitivity_table` with
`\tableref{cash_transfers_sensitivity_table}`, and to
`OHIE-graphics` with `\figref{OHIE-graphics}`.  Let's not use `ref`.
