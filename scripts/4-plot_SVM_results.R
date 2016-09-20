# Jaclyn Taroni Sep 2016
# Purpose: Plot SVM results
# Usage: Rscript scripts/4-plot_SVM_results.R

#### packages & source ---------------------------------------------------------

library(dplyr)
library(reshape2)
library(ggplot2)

source("util/theme_black.R")

#### setup files ---------------------------------------------------------------

# directories
res.dir <- "results"
plot.dir <- "plots"
data.dir <- "data"

# input
svm.file <- 
  paste0(res.dir, "/tubulointerstitium_pos_glomeruli_neg_kidney_SVM_output.tsv")
tub.file <- paste0(res.dir, "/tubulointerstitium_DEGS_LN_v_LD.tsv")
glom.file <- paste0(res.dir, "/glomeruli_DEGS_LN_v_LD.tsv")
symbol.entrez.file <- paste0(data.dir, "/symbol_entrez.txt")

# output
scatter.file <- paste0(plot.dir, "/t-statistic_SVM.score_scatterplot.png")
symbol.svm.file <- 
  paste0(res.dir, 
  "/tubulointerstitium_pos_glomeruli_neg_kidney_SVM_output_gene_symbols.tsv")

#### read in data --------------------------------------------------------------

tub.degs <- read.delim(tub.file)
glom.degs <- read.delim(glom.file)
svm.df <- read.delim(svm.file, header = F)
colnames(svm.df) <- c("Entrez.ID", "Standard", "SVM.Score")
symbol2entrez.mapping <- read.delim(symbol.entrez.file, header = F)
colnames(symbol2entrez.mapping) <- c("Gene.Symbol", "Entrez.ID")

#### map gene identifiers back to gene symbols ---------------------------------

svm.symbol.df <- dplyr::right_join(symbol2entrez.mapping,
                                   svm.df,
                                   by = "Entrez.ID",
                                   copy = TRUE)

write.table(svm.symbol.df, symbol.svm.file,
            quote = F, row.names = F, col.names = F, sep = "\t")

#### scatterplot ---------------------------------------------------------------

mstr.df <- dplyr::inner_join(dplyr::select(tub.degs, Gene.Symbol, t),
                             dplyr::select(glom.degs, Gene.Symbol, t),
                             by = "Gene.Symbol")
mstr.df <- dplyr::inner_join(mstr.df,
                             dplyr::select(svm.symbol.df, Gene.Symbol, Standard,
                                           SVM.Score),
                             by = "Gene.Symbol")
colnames(mstr.df) <- c("Gene.Symbol", "Tub.t", "Glom.t", 
                       "Standard", "SVM.Score")
mstr.df$Standard <- sub("1", "Tub", mstr.df$Standard)
mstr.df$Standard <- sub("-Tub", "Glom", mstr.df$Standard)
mstr.df$Standard <- sub("0", "None", mstr.df$Standard)
mstr.df$Standard <- as.factor(mstr.df$Standard)

png(scatter.file, units = "in", width = 12, height = 12, res = 300)
ggplot2::ggplot(mstr.df, aes(x = Tub.t, y = Glom.t, 
                             colour = SVM.Score, shape = Standard)) +
  geom_point(alpha = 1/5) +
  scale_color_gradient2(low = "blue", high = "yellow") +
  theme_black() +
  ggtitle("Lupus Nephritis Kidney") +
  labs(x = "tubulointerstitium t-statistic",
       y = "glomeruli t-statistic") +
  theme(legend.key = element_rect(fill = "white"))
dev.off()
