# Lupus-NetSVM
Use machine learning on GIANT tissue-specific functional genomic networks to 'boost' disease differential expression signal.

## Overview

The objective of this project is to identify changes in lupus nephritis kidney
that are *highly-specific* to the tubulointerstitium and glomeruli. 

We used gene expression from lupus nephritis kidney biopsies originally analyzed 
as part of [Berthier et al. J Immunol. 2012.](https://www.ncbi.nlm.nih.gov/pubmed/22723521). 
This expression data is deposited in GEO under the accession number
[GSE32591](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32591).

We identify genes that are upregulated in lupus nephritis kidneys 
(differentially expressed genes, abbreviated DEGS) compared to living donor 
controls in each compartment separately. We then supply these genes as negative 
(glomeruli DEGS) or positive (tubulointerstitium DEGS) standards to a linear
SVM that learns the genes' connectivity patterns in the kidney functional 
genomic network from [GIANT](giant.princeton.edu) 
([Greene et al. Nat Genet. 2015](http://www.ncbi.nlm.nih.gov/pubmed/25915600)).
