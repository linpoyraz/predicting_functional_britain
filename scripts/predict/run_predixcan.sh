#!/bin/bash

module load python/2

tiss=$1
dir=$2
out=$3

mkdir -p ${out}/predictions/${tiss}

/project/mathilab/ppoyraz/time_series_predixcan/bin/PrediXcan.py --predict --dosages ${dir} --dosages_prefix chr --samples sample.txt --weights /project/mathilab/colbranl/data/zhou2020_JTI_models/${tiss}.db --output_dir ${out}/predictions/${tiss}/
