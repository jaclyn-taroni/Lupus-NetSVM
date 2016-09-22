#!/bin/sh

scripts/0-data_download.sh
Rscript scripts/1-get_differentially_expressed_genes.R
Rscript scripts/2-generate_standards.R
scripts/3-run_SVM.sh
Rscript scripts/4-plot_SVM_results.R
