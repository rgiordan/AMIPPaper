This directory contains scripts for running and analyzing the microcredit
ADVI example.

Unfortunately all ADVI analysis, including using `python` and `R`, is contained
in the `python` folder.  I don't want to rename it now because I don't want to
break anything.  Let's rename it when we do the final run.

Fitting in Python:
- `microcredit_vb.ipynb`: Run and save an initial fit.
- `sensitivity.ipynb`: Load the initial fit, compute and save sensitivity metrics.
- `refit.py`, `run_refits.py`: Re-fit the model with new weights to check the
sensitivity results.

Analysis in R:
- `process_results_for_paper.R`: Generate a single, relatively small `Rdata`
file summarizing the original fit.
- `generate_refit_weights.R`:  For select parameters, save weights in a
format that can be used by `refit.py`.
- `generate_graphs.R`: An example of loading the results and creating the sort
of graphs that might go in the paper.
- `prior_sensitivity.R`: Investigate prior sensitivity.

Other, less important scripts:
- `checking.R`: Systematically test that the Stan and Python models are
identical.
- `analyze_smuber.R`: Dig into what the Smuber loss is doing.
