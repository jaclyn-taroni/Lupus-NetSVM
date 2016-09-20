#!/bin/sh

SVMperfer -l results/tubulointerstitium_pos_glomeruli_neg_DEGS_standard.tsv -o results/tubulointerstitium_pos_glomeruli_neg_kidney_SVM_output.tsv -i data/kidney.dab -a -k 0 -c 5 -e 2 -t 50

