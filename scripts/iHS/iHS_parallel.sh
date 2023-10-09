#!/bin/bash

# Lin Poyraz 
# This script runs iHS generation in parallel 

module load R

results=/project/mathilab/ppoyraz/time_series_predixcan/results/data/scan/iHS_split 
regression=./output/samples01_time_series_regression_results.txt
final=/project/mathilab/ppoyraz/time_series_predixcan/results/data/scan/iHS_results.txt

rm ${final}
touch ${final}

mkdir -p ${results}
mkdir -p temp

for chr in {1..22}; do
	bsub "zgrep ^chr${chr}$'\t' ${regression} > ./temp/${chr} 	
	Rscript iHS_results.R ./temp/${chr} 
	rm ./temp/${chr}"
done

cat ${results}/1_results.txt > ${final} 

for chr in {2..22};do
	tail -n +2 ${results}/${chr}_results.txt >> ${final}
done 



