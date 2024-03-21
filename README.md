# AdversarialInfluenceWorkbench

This repository contains code to reproduce our paper,
[An Automatic Finite-Sample Robustness Metric: When Can Dropping a Little Data Make a Big Difference?](https://arxiv.org/abs/2011.14999) by Tamara Broderick, Ryan Giordano, and Rachael Meager.

The writing directory is `writing/output/`.  There is data processing
code in both `examples` and in `writing/applications`.

The best guide to reproducing the paper and its analyses is found in
`writing/output/makefile`.  To run it, first

1) Set the `GIT_REPO_LOC` variable at the top of `makefile` to point to the full path of the
   location of the cloned `AMIPPaper` repository
2) Run `make all` in the `output` directory.
3) Follow the instructions to download the needed data.
4) Continue to run `make all` and follow the instructions until the paper succesfully compiles.

To clear paper output, run `make clean`.

To pre-process data for individual analyses, you can run any of the subsidiary targets:

- `make sim_data`
- `make cash_data`
- `make mc_data`
- `make mc_data`

Note that the most complicated analysis, the mixture model, requires more
effort to run than the `R` analysis.

In order to simply compile the paper without running any of the individual
analyses, you can uncompress the file `writing/output/applications_data.tgz`,
which contains the processed data we used for our original paper.
The contents of the arxiv should replace the contents of the
`writing/output/applications_data/` directory, and satisfy
the `makefile` target `$(PP_DATA)`.
These can also be used to sanity check your own runs against ours.


If you have problems reproducing any aspect of the pipeline, please
send Ryan an email.