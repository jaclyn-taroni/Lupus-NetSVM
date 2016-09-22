# Lupus-NetSVM

**Use machine learning on GIANT tissue-specific functional genomic networks to 
'boost' disease differential expression signal.**

A similar approach has been applied to drug studies with molecular data in
the autoimmune disease systemic sclerosis; that work will be presented as an 
oral presentation at the *2016 American College of Rheumatology Annual Meeting*.

## Rationale

We aim to put differentially expressed genes into biological context in a manner 
that does not rely on curated gene sets. This allows us to identify relationships 
that are unknown put supported by data and to do so in a manner that is 
context (tissue) specific. This approach has been succesfully applied to summary 
statistics in GWAS data in a method termed **NetWAS**. See 
[Greene et al. Nat Genet. 2015](http://www.ncbi.nlm.nih.gov/pubmed/25915600).

## Overview

The objective of this project is to identify changes in lupus nephritis kidney
that are *highly-specific* to the tubulointerstitium and glomeruli. 

We used gene expression from lupus nephritis kidney biopsies originally analyzed 
as part of 
[Berthier et al. J Immunol. 2012](https://www.ncbi.nlm.nih.gov/pubmed/22723521). 
This expression data is deposited in GEO under the accession number
[GSE32591](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32591).

We identify genes that are upregulated in lupus nephritis kidneys 
(differentially expressed genes, abbreviated DEGS) compared to living donor 
controls in each compartment separately. We then supply these genes as negative 
(glomeruli DEGS) or positive (tubulointerstitium DEGS) standards to a linear
SVM that learns the genes' connectivity patterns in the kidney functional 
genomic network from [GIANT](giant.princeton.edu) 
([Greene et al. Nat Genet. 2015](http://www.ncbi.nlm.nih.gov/pubmed/25915600)).

## Analysis

Scripts for our analyses are in the scripts folder. 
The entire pipeline can be performed by running
```
./run_experiments.sh
```
in the top directory.


## Dependencies

R, the R packages in `install_packages.R` and 
[Sleipnir](http://libsleipnir.bitbucket.org/index.html) with SVMPerf.

We are in the process of preparing a Docker image. 
