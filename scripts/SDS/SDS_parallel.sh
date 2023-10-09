#!/bin/bash

# this script runs SDS in parallel 
results=/project/mathilab/ppoyraz/time_series_predixcan/results/data/scan/SDS_split 
regression=/project/mathilab/ppoyraz/time_series_predixcan/results/data/expression/time_series_regression_results.txt
final=/project/mathilab/ppoyraz/time_series_predixcan/results/data/scan/SDS_results.txt


mkdir -p ${results}
mkdir -p temp

for chr in {1..22}; do
	bsub "zgrep ^chr${chr}$'\t' ${regression} > ./temp/${chr} 
	Rscript SDS_results.R ./temp/${chr} 
	rm ./temp/${chr}"
done 

cat ${results}/1_results.txt > ${final} 

for chr in {2..22};do
	tail -n +2 ${results}/${chr}_results.txt >> ${final}
done 


