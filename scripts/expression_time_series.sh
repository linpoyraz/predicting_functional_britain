#!/bin/bash

# This script runs transcriptome-wide selection scan 

# make a directory for the output 
mkdir -p output 

#Rscript p_cutoff_time_series_regression.R ./data ./output samples_01

Rscript random_time_series_regression.R ./data


