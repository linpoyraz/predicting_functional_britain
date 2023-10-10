#!/bin/bash

# Lin Poyraz 
# This script generates mean imputation scores for all genes. 

results=./output/imputation_results_split
regression=./output/samples01_time_series_regression_results.txt
final=./output/imputation_results.txt

mkdir -p ${results}
mkdir -p temp

for chr in {1..22}; do
	zgrep ^chr${chr}$'\t' ${regression} > ./temp/${chr} 
	Rscript imputation_results.R ./temp/${chr} 
	rm ./temp/${chr}
done 

rm ${final} 

cat ${results}/1_results.txt > ${final} 

for chr in {2..22};do
	tail -n +2 ${results}/${chr}_results.txt >> ${final}
done 


