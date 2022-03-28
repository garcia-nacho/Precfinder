#!/bin/bash

source activate pangolin 
pangolin --update
pangolin --usher -t 10 ${1} --outfile ${1}_pango.csv
conda deactivate 

