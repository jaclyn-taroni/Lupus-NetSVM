#!/bin/sh

data_folder='data/'
wget -i data/expression_data_urls.txt '--directory-prefix='$data_folder
cd data
gunzip *.gz
wget http://giant.princeton.edu/static/networks_dab/kidney.dab
