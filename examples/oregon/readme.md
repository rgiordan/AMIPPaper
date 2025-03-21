# README

prepare_OHIE_data.R is an R script that replicates the authors' .do file. If write_to_disk is set to true then a csv is produced named data_for_analysis_R.csv 


QJE_Tables.R replicates tables 5, 6, 8, 9 and 10. The initial tables I didn't bother replicating as they're either summary stats or first stages. Table 11 should replicate but I haven't got around to writing the code (it's not a very interesting table either, imo). The remaining tables draw on data not publicly available. The tables that do replicate produce perfect point estimates and standard errors.

I use R 3.5.3 and packages used can be found in install.R