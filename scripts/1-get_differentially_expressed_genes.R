# Jaclyn Taroni Sep 2016
# Some code adapted from GEO2R R script generated Sep 19 2016
# Purpose: Identify genes that are differentially expressed between control
# and lupus nephritis samples in two compartments of the kidney -- 
# glomeruli and tubulointerstitium
# Usage: Rscript scripts/1-get_differentially_expressed_genes.R

#### packages ------------------------------------------------------------------

suppressMessages(library(checkpoint))
suppressMessages(checkpoint("2016-09-19", checkpointLocation = "."))

library(Biobase)
library(limma)
library(data.table)
library(dplyr)

#### custom functions ----------------------------------------------------------

BuildExpressionSet <- function(exprs.df) {
  # This function takes a data.frame of expression data and builds an
  # ExpressionSet from it.
  # 
  # Args:
  #   exprs.df: a data.frame of expression data, rows are genes, columns are 
  #             samples; already processed to remove affy internal controls,
  #             genes below baseline expression (and the corresponding column)
  # 
  # Returns: 
  #   eset: an ExpressionSet with duplicate genes removed
  #             
  
  exprs.mat <- as.matrix(exprs.df[, 4:ncol(exprs.df)])
  rownames(exprs.mat) <- exprs.df$`Gene Symbol`
  # for duplicate genes, remove second instance
  to.rm <- which(duplicated(rownames(exprs.mat)))
  exprs.mat <- exprs.mat[-to.rm, ]
  eset <- Biobase::ExpressionSet(assayData = exprs.mat)
  return(eset)
  
}

GetDEGS <- function(gset, grp.vector) {
  # This function takes an ExpressionSet and a vector that indicates the two
  # groups to be compared for differential expression.
  # Adapted from GEO2R R Script generated Sep 19 2016.
  # 
  # Args:
  #   gset: an ExpressionSet
  #   grp.vector: a vector of group labels
  #   
  # Returns:
  #   tT: a data.frame that contains the t-statistic, p-value, FDR-adjusted
  #       p-value and gene symbol of all genes tested
  #     
  
  # Error-handling
  if (dim(gset)[[2]] != length(grp.vector)) {
    stop("grp.vector length and the number of samples in gset must be equal")
  }
  if (class(grp.vector) != "factor") {
    grp.vector <- as.factor(grp.vector)
  }
  if (!(all(levels(grp.vector) == c("G0", "G1")))) {
    stop("grp.vector levels must be G0 and G1")
  }
  
  gset$description <- grp.vector
  design <- model.matrix(~ description + 0, gset)
  colnames(design) <- levels(grp.vector)
  fit <- limma::lmFit(gset, design)
  cont.matrix <- limma::makeContrasts(G1-G0, levels = design)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, 0.01)
  tT <- limma::topTable(fit2, adjust = "fdr", sort.by = "t", 
                 number = dim(gset)[[1]])
  tT <- subset(tT, select=c("t","adj.P.Val","P.Value"))
  tT <- dplyr::bind_cols(data.frame(rownames(tT)), tT)
  colnames(tT)[1] <- "Gene.Symbol"
  
  return(tT)
  
}

#### file setup ----------------------------------------------------------------

# directories
data.dir <- "data"
res.dir <- "results"

# input files
tub.file <- paste0(data.dir, "/GSE32591_matrix_tubulointerstitium.txt")
glom.file <- paste0(data.dir, "/GSE32591_matrix_glomeruli.txt")

# degs output files 
tub.degs.file <- paste0(res.dir, "/tubulointerstitium_DEGS_LN_v_LD.tsv")
glom.degs.file <- paste0(res.dir, "/glomeruli_DEGS_LN_v_LD.tsv")

#### read in data & filter -----------------------------------------------------

tub.df <- data.table::fread(tub.file, data.table = FALSE)
glom.df <- data.table::fread(glom.file, data.table = FALSE)

# remove genes below baseline and Affymetrix Internal Controls
tub.internal <- grep("Affymetrix Internal control", tub.df$Description)
glom.internal <- grep("Affymetrix Internal control", glom.df$Description)
tub.below.baseline <- 
  which(tub.df$`Expressed above the defined expression baseline` == "No")
glom.below.baseline <- 
  which(glom.df$`Expressed above the defined expression baseline` == "No")
tub.df <- tub.df[-union(tub.below.baseline, tub.internal), 
                   -c(4, as.integer(which(apply(tub.df, 2, 
                                                function(x) all(is.na(x))))))]
glom.df <- glom.df[-union(glom.below.baseline, glom.internal), -4]

rm(tub.below.baseline, glom.below.baseline)


# build ExpressionSets
tub.eset <- BuildExpressionSet(tub.df)
glom.eset <- BuildExpressionSet(glom.df)
rm(tub.df, glom.df)

#### differential expression analysis ------------------------------------------

# classes / group setup
tub.grp <- rep("0", dim(tub.eset)[[2]])
tub.grp[grep("LD", colnames(exprs(tub.eset)))] <- "1"
tub.grp <- as.factor(paste0("G", tub.grp))

glom.grp <- rep("0", dim(glom.eset)[[2]])
glom.grp[grep("LD", colnames(exprs(glom.eset)))] <- "1"
glom.grp <- as.factor(paste0("G", glom.grp))

# get tables of differentially expressed genes
tub.degs <- GetDEGS(gset = tub.eset, grp.vector = tub.grp)
glom.degs <- GetDEGS(gset = glom.eset, grp.vector = glom.grp)

# write DEGS to file
write.table(tub.degs, file = tub.degs.file, sep = "\t", row.names = F, 
            quote = F)
write.table(glom.degs, file = glom.degs.file, sep = "\t", row.names = F, 
            quote = F)
