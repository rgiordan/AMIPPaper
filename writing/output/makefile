# Set this variable to the full path of your copy of the paper's git repo
GIT_REPO_LOC := /home/rgiordan/Documents/git_repos/AMIPPaper



###################################
# You should not have to edit anything below this line.

# Files for the paper
PAPER_DIR := $(GIT_REPO_LOC)/writing/output

# Postprocessed data ready to go in the paper
DATA_DIR := $(PAPER_DIR)/applications_data

# Sensitivity analyses
APP_DIR := $(GIT_REPO_LOC)/writing/applications

# The original analyses
EXAMPLES_DIR := $(GIT_REPO_LOC)/examples

# The name of our conda environment
CONDA_ENV := amip-paper-2023


###################################
# Set up a conda environment.

.PHONY: conda_env
conda_env:
	conda env create -f $(EXAMPLES_DIR)/microcredit_mixture/environment.yml

# R packages
.PHONY: r_packages
r_packages:
	Rscript $(PAPER_DIR)/figures_knitr/install.R



###################################
# Data

SIM_DATA := $(DATA_DIR)/simulations/noise_grid.Rdata
CASH_DATA := $(DATA_DIR)/cash_transfers/cash_transfers_results.Rdata
OHIE_DATA := $(DATA_DIR)/ohie/OHIE_results.Rdata
MC_DATA := $(DATA_DIR)/microcredit/microcredit_results.Rdata
MCT_DATA := $(DATA_DIR)/microcredit/microcredit_temptation_results.Rdata
MCMIX_DATA := $(DATA_DIR)/microcredit_mixture/microcredit_mixture_results.Rdata




# Simulations
$(SIM_DATA):
	Rscript $(APP_DIR)/simulations/ols_example.R --git_repo_loc=$(GIT_REPO_LOC)
.PHONY: sim_data
sim_data: $(SIM_DATA)	

.PHONY: clean_sim_data
clean_sim_data:
	rm -f $(SIM_DATA)


# Cash transfers
CASH_DTA := $(EXAMPLES_DIR)/examples/cash_transfers/angelucci-digiorgi-data-code/table1.dta
$(CASH_DTA):
	echo The cash transfers analysis requires $(CASH_DTA)
	echo The cash transfers raw data can be found here:
	echo https://www.openicpsr.org/openicpsr/project/113289/version/V1/view
	exit 1

CASH_RDS := $(EXAMPLES_DIR)/cash_transfers/processed_CT_data.rds
$(CASH_RDS):
	Rscript $(APP_DIR)/cash_transfers/read-cash-data-into-R.R --git_repo_loc=$(GIT_REPO_LOC)

$(CASH_DATA): $(CASH_RDS)
	Rscript $(APP_DIR)/cash_transfers/cash_consumption_sensitivity.R --git_repo_loc=$(GIT_REPO_LOC)
.PHONY: cash_data
cash_data: $(CASH_DATA)

.PHONY: clean_cash_data
clean_cash_data:
	rm -f $(CASH_RDS)
	rm -f $(CASH_DATA)


# OHIE
OHID_DATA_DIR := $(EXAMPLES_DIR)/oregon/Data/OHI/OHIE_Public_Use_Files/OHIE_Data
OHIE_CSV := $(OHID_DATA_DIR)/data_for_analysis_R.csv

# Downloaded data
$(OHID_DATA_DIR)/oregonhie_descriptive_vars.dta:
	echo The Oregon medicaid analysis requires \
		$(OHID_DATA_DIR)/oregonhie_descriptive_vars.dta
	echo All Oregon medicaid data is from here:
	echo https://www.nber.org/research/data/oregon-health-insurance-experiment-data
	exit 1
$(OHID_DATA_DIR)/oregonhie_stateprograms_vars.dta:
	echo The Oregon medicaid analysis requires \
		$(OHID_DATA_DIR)/oregonhie_stateprograms_vars.dta
	echo All Oregon medicaid data is from here:
	echo https://www.nber.org/research/data/oregon-health-insurance-experiment-data
	exit 1

$(OHIE_CSV): $(OHID_DATA_DIR)/oregonhie_descriptive_vars.dta $(OHID_DATA_DIR)/oregonhie_stateprograms_vars.dta
	Rscript $(APP_DIR)/OHIE/prepare_OHIE_data.R --git_repo_loc=$(GIT_REPO_LOC)
$(OHIE_DATA): $(OHIE_CSV)
	Rscript $(APP_DIR)/OHIE/OHIE-sensitivity-table-9-with-IV.R --git_repo_loc=$(GIT_REPO_LOC)
.PHONY: ohie_data
ohie_data: $(OHIE_DATA)

.PHONY: clean_ohie_data
clean_ohie_data:
	rm -f $(OHIE_CSV)
	rm -f $(OHIE_DATA)



# Microcredit

MC_BASE_RDATA := $(EXAMPLES_DIR)/microcredit/microcredit_project_data.RData
$(MC_BASE_RDATA):
	The microcredit analysis requires $(MC_BASE_RDATA)
	The microcredit project repo officially is here:
	echo https://www.openicpsr.org/openicpsr/project/155821/version/V1/view 
	echo The data and associated data construction script (for documentation purposes) is here:
	echo https://www.openicpsr.org/openicpsr/project/155821/version/V1/view?path=/openicpsr/155821/fcr:versions/V1/aggregating-dist-TEs/data&type=folder
	exit 1

# Something is messed up here
MC_RDATA := $(APP_DIR)/microcredit/microcredit_profit_data.Rdata
MCT_RDATA := $(APP_DIR)/microcredit/microcredit_temptation_data.Rdata

$(MC_RDATA):
	Rscript $(APP_DIR)/microcredit/preprocess_microcredit_data.R --git_repo_loc=$(GIT_REPO_LOC)
$(MCT_RDATA):
	Rscript $(APP_DIR)/microcredit/preprocess_microcredit_data.R --git_repo_loc=$(GIT_REPO_LOC)

$(MC_DATA): $(MC_RDATA)
	Rscript $(APP_DIR)/microcredit/microcredit_profit_sensitivity.R --git_repo_loc=$(GIT_REPO_LOC)

$(MCT_DATA): $(MCT_RDATA)
	Rscript $(APP_DIR)/microcredit/microcredit_temptation_sensitivity.R --git_repo_loc=$(GIT_REPO_LOC)
.PHONY: mc_data
mc_data: $(MC_DATA) $(MCT_DATA)

.PHONY: clean_mc_data
clean_mc_data:
	rm -f $(MC_RDATA)
	rm -f $(MCT_RDATA)
	rm -f $(MC_DATA)
	rm -f $(MCT_DATA)



# Microcredit mixture

# Note: this is still a work in progress.

# Process the output of the notebooks
# microcredit_vb.ipynb and sensitivity.ipynb, and the R script run_refits.R
# into a format that we can use in the paper.

SHELL = /bin/bash
CONDA_ACTIVATE := eval "$$(conda shell.bash hook)"; source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

MCMIX_DRAWS := 30

MCMIX_DIR := $(EXAMPLES_DIR)/microcredit_mixture
MCMIX_CSV := $(MCMIX_DIR)/R/data/microcredit_project_data_cleaned.csv

# This generates both the csv and the MCMC samples unfortunately.
$(MCMIX_CSV): $(MC_BASE_RDATA)
	$(MCMIX_DIR)/R/preprocess_microcredit_data_for_stan.R \
		--git_repo_loc=$(GIT_REPO_LOC)

# MCMC samples:
# output/microcredit_profit_lognormal_tailored_hierarchical_pdf_output_4000_iters.RData?
# microcredit_profit_lognormal_tailored_hierarchical_pdf_output_5000_iters.RData?
# MCMIX_MCMC := $(EXAMPLE_DIR)/microcredit_mixture/data/...

.ONESHELL:
MCMIX_INITIAL_FIT := $(MCMIX_DIR)/python/output/microcredit_project_advi_$(MCMIX_DRAWS)draws_smuberFalse.npz
$(MCMIX_INITIAL_FIT): $(MCMIX_CSV)
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	$(MCMIX_DIR)/python/microcredit_vb.py --git_repo_loc=$(GIT_REPO_LOC) --num_draws=$(MCMIX_DRAWS)
.PHONY: mcmix_initial_fit
mcmix_initial_fit: $(MCMIX_INITIAL_FIT)


.ONESHELL:
MCMIX_SENS := $(MCMIX_DIR)/python/output/microcredit_project_weight_sensitivity_$(MCMIX_DRAWS)draws_smuberFalse.npz
$(MCMIX_SENS): $(MCMIX_INITIAL_FIT) $(MCMIX_CSV)
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	$(MCMIX_DIR)/python/sensitivity.py --git_repo_loc=$(GIT_REPO_LOC) --num_draws=$(MCMIX_DRAWS)
.PHONY: mcmix_sens
mcmix_sens: $(MCMIX_SENS)


# Requires microcredit_vb_%s_combined_results.Rdata and microcredit_vb_%s_refit_results.Rdata
MCMIX_SUMMARY_DATA := $(MCMIX_DIR)/data/microcredit_vb_draws$(MCMIX_DRAWS)_smuberFalse_combined_results.Rdata
.ONESHELL:
$(MCMIX_SUMMARY_DATA): $(MCMIX_INITIAL_FIT) $(MCMIX_SENS) $(MCMIX_CSV)
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	Rscript $(MCMIX_DIR)/python/save_initial_fit_as_rdata.R \
		--git_repo_loc=$(GIT_REPO_LOC) \
		--num_draws=$(MCMIX_DRAWS) \
		--python_path=$$(which python)


MCMIX_REFIT_DATA := $(MCMIX_DIR)/data/microcredit_vb_draws$(MCMIX_DRAWS)_smuberFalse_refit_results.Rdata
.ONESHELL:
$(MCMIX_REFIT_DATA): $(MCMIX_SUMMARY_DATA) $(MCMIX_INITIAL_FIT) $(MCMIX_SENS) $(MCMIX_CSV)
	echo Running refits.  This will probably take a while.
	echo However, by default, terminating and re-running this command will start off where it left off.
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	Rscript $(MCMIX_DIR)/python/run_refits.R \
		--git_repo_loc=$(GIT_REPO_LOC) \
		--num_draws=$(MCMIX_DRAWS) \
		--python_path=$$(which python)
.PHONY: mcmix_refits
mcmix_refits: $(MCMIX_REFIT_DATA) 


.ONESHELL:
$(MCMIX_DATA): $(MCMIX_SUMMARY_DATA) $(MCMIX_REFIT_DATA)
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	Rscript $(APP_DIR)/microcredit_mixture/process_results_for_paper.R \
		--git_repo_loc=$(GIT_REPO_LOC) \
		--num_draws=$(MCMIX_DRAWS)
.PHONY: mcmix_data
mcmix_data: $(MCMIX_DATA) 


.PHONY: clean_mcmix_initial_fit
clean_mcmix_initial_fit:
	rm -f $(MCMIX_INITIAL_FIT)

.PHONY: clean_mcmix_refits
clean_mcmix_refits:
	rm -f $(MCMIX_DIR)/python/reweightings/amis_weight_draws$(MCMIX_DRAWS)_*smuberFalse.npz
	rm -f $(MCMIX_DIR)/python/reweightings/refit_*_$(MCMIX_DRAWS)draws.npz

.PHONY: clean_expensive_mcmix_data
clean_expensive_mcmix_data: clean_mcmix_refits clean_mcmix_initial_fit


.PHONY: clean_mcmix_data
clean_mcmix_data:
	rm -f $(MCMIX_DATA)
	rm -f $(MCMIX_SUMMARY_DATA)
	rm -f $(MCMIX_SENS)
	rm -f $(MCMIX_REFIT_DATA)
	echo Not removing very time-consuming MCMIX analyses \(initial fit and refits\).
	echo To remove the initial MCMIX fit and refits, run 'make clean_expensive_mcmix_data'


#~/Documents/git_repos/AdversarialInfluenceWorkbench/examples/
.PHONY: clean_data
clean_data: clean_sim_data clean_cash_data clean_ohie_data clean_mc_data clean_mcmix_data


###################################
# The actual paper.

.PHONY: all
all: main.pdf

main.pdf: figures_knitr.tex
	pdflatex main.tex
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex

PP_DATA := $(SIM_DATA) $(CASH_DATA) $(OHIE_DATA) $(MC_DATA) $(MCT_DATA) $(MCMIX_DATA)

figures_knitr.tex: $(PP_DATA)
		Rscript -e 'library(knitr); knit("figures_knitr.Rnw")'


.PHONY: recompile_knitr
recompile_knitr:
		Rscript -e 'library(knitr); knit("figures_knitr.Rnw")'
		pdflatex main.tex


###################################
# Cleaning



.PHONY: clean_tex
clean_tex:
		for latextext in aux bbl blg log out pdf ; do \
				rm -f main.$$latextext; \
		done

.PHONY: clean_knitr
clean_knitr:
		rm -f figures_knitr.tex
		rm -f figure/*
		rm -f cache/*

.PHONY: clean
clean: clean_knitr clean_tex