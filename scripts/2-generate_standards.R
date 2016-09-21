# Jaclyn Taroni Sep 2016
# Purpose: Generate standards for SVMperfer; we are interested in learning more
# about lupus nephritis changes specific to each compartment
# Usage: Rscript scripts/2-generate_standards.R

#### packages ------------------------------------------------------------------

suppressMessages(library(checkpoint))
suppressMessages(checkpoint("2016-09-19", checkpointLocation = "."))

library(dplyr)

#### setup filenames -----------------------------------------------------------
# directories 
res.dir <- "results"
data.dir <- "data"

# input
tub.file <- paste0(res.dir, "/tubulointerstitium_DEGS_LN_v_LD.tsv")
glom.file <- paste0(res.dir, "/glomeruli_DEGS_LN_v_LD.tsv")
symbol.entrez.file <- paste0(data.dir, "/symbol_entrez.txt")

# output: standard
standard.file <- 
  paste0(res.dir, "/tubulointerstitium_pos_glomeruli_neg_DEGS_standard.tsv")

#### read in data --------------------------------------------------------------

tub.degs <- read.delim(tub.file)
glom.degs <- read.delim(glom.file)
symbol2entrez.mapping <- read.delim(symbol.entrez.file, header = F)
colnames(symbol2entrez.mapping) <- c("Gene.Symbol", "Entrez.ID")

#### generate standard --------------------------------------------------------

# what genes are higher in lupus nephritis samples?
tub.LN.genes <- 
  as.character(tub.degs$G[which(tub.degs$P < 0.001 & tub.degs$t < 0)])
glom.LN.genes <-
    as.character(glom.degs$G[which(glom.degs$P < 0.001 & glom.degs$t < 0)])

# remove the genes that are differentially expressed in both compartments
overlap.genes <- intersect(tub.LN.genes, glom.LN.genes)
tub.LN.genes <- setdiff(tub.LN.genes, overlap.genes)
glom.LN.genes <- setdiff(glom.LN.genes, overlap.genes)

# build standard -- want to learn differences between lupus nephritis in
# the different compartments
standard.df <- cbind(c(tub.LN.genes, glom.LN.genes),
                     c(rep(1, length(tub.LN.genes)), 
                       rep(-1, length(glom.LN.genes))))
colnames(standard.df) <- c("Gene.Symbol", "Standard")
standard.df <- dplyr::right_join(symbol2entrez.mapping,
                                 standard.df,
                                 by = "Gene.Symbol",
                                 copy = TRUE)
# remove genes that do not map to Entrez ID and gene symbol identifier column
standard.df <- standard.df[complete.cases(standard.df), -1]

# output standard file 
write.table(standard.df, standard.file, quote = F, row.names = F, col.names = F,
            sep = "\t")
