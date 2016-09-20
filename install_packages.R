# Jaclyn Taroni Sep 2016
# Usage: This script is to be run upon docker image build

library("methods")

mirror <- "http://cran.us.r-project.org"
install.packages('checkpoint', repos = mirror)

library("checkpoint")

#### install cran packages -----------------------------------------------------

dir.create('.checkpoint')
checkpoint('2016-09-19', checkpointLocation = '.')
cran.pkgs <- c("data.table", "ggplot2", "dplyr", "reshape2", "scales", 
                "gridExtra")
install.packages(cran.pkgs)

#### install cran packages -----------------------------------------------------

source("https://bioconductor.org/biocLite.R")
bioconductor.pkgs <- c("Biobase","limma")
biocLite(bioconductor.pkgs, suppressUpdates = TRUE)
