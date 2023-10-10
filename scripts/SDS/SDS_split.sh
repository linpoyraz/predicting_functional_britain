#!/bin/bash

# Lin Poyraz 
# This script splits Field et al. 2016 UK10K SDS results by chromosome

my_dir=/project/mathilab/ppoyraz/scan/SDS
mkdir -p ${my_dir}

for chr in {1..22};do
  echo "Preparing chromosome ${chr} files"
  zgrep ^${chr}$'\t' /project/mathilab/ppoyraz/scan/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz |  sort -nk2 > ${my_dir}/SDS_${chr}.tab
  gzip ${my_dir}/SDS_${chr}.tab
done 
