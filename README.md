# Predicting functional consequences of recent natural selection in Britain 

## figures/

This directory includes all files and scripts required to generate all main figures in the manuscript. 

- `main_figures.Rmd` and `main_figures.html` \
This markdown generates all main figures, as well as Supplementary Figure 3 in the manuscript. 

### figures/data/

This directory includes all output files required to generate all main figures in the manuscript. 

- `samples01_time_series_regression_results.txt.gz` \
The results of the transcriptome-wide selection scan. This file is produced by the script `expression_time_series.sh.`

- `samples01_random_time_series_regression_results.txt.gz` \
The results of the transcriptome-wide selection scan with date labels randomized. This file is produced by the script `expression_time_series.sh.`

- `s_scan_all_brit2_01_window.txt.gz` \
The output of the window-based genome-wide selection scan as described by [Mathieson and Terhorst 2022](https://genome.cshlp.org/content/32/11-12/2057). 

- `s_scan_all_brit2_01_random_window.txt.gz` \
The output of the window-based genome-wide selection scan with randomized date labels as described by [Mathieson and Terhorst 2022](https://genome.cshlp.org/content/32/11-12/2057). 

- `model_info.txt.gz` \
Tissue information for the JTI models with the highest training $R^2$. 

- `imputation_results.txt.gz` \
The mean imputation qualities of each gene. This file is the output of the script `imputation_results.sh.` 

- `iHS_results.txt.gz` \
Gene-level iHS scores. This file is the output of `iHS_parallel.sh`.

- `SDS_results.txt.gz` \
Gene-level SDS scores. This file is the output of `SDS_parallel.sh`.

### figures/data/allele_frequency/

This directory includes information regarding all SNPs included in the JTI model for all FDR significant transcriptome-wide scan genes. These files are the outputs of `allele_frequencies.R`. 

### figures/data/expression/

This directory includes the predicted expression of each FDR significant transcriptome-wide scan gene for all samples used in this manuscript. These files are the outputs of `expression.R`

## scripts/

- `allele_frequencies.R` \
This script outputs dosage information and JTI weights for all SNPs included in the JTI model for the specified gene. Paths are absolute with regard to the cluster, so it won't be possible to run this script on another machine without significant editing.

- `expression.R` \
This script outputs predicted expression of each FDR significant transcriptome-wide gene for all samples used in the manuscript. 

### scripts/predict/

This directory includes all scripts required to generate predicted gene expression values for all samples included in this manuscript. Paths are absolute with regard to the cluster, so it won't be possible to run this script on another machine without significant editing.

- `create_dosages.sh` \
This script prepares ancient and 1K Genomes GBR files. It calls `vcf2dosageprob.py` to convert vcf files to dosage files, and `update_rsIDs.py` which updates rsIDs in a dosage file to the ones in PrediXcan model by genomic location. 

- `prepare_predixcan_dir.R` \
This script adds sample.txt file to the Predixcan directory. 

- `run_predixcan.sh` \
This script calls `PrediXcan.py` and generates predicted gene expression values for all samples. 

### scripts/iHS/

This directory includes scripts required to generate gene-level iHS scores. Paths are absolute with regard to the cluster, so it won't be possible to run these scripts on another machine without significant editing. The output of these scripts is included in `figures/data/iHS_results.txt.gz`. 

- `iHS_generate.sh` \
This script polarizes 1kg GBR variants using PanTro6 using `ancestral_match.R`, calls `interpolate_map.R` to interpolate recombination maps, and generates iHS scores using selscan. It then calls `iHS_split.R` to add ref/alt and rsid information to the normalized output. 

- `iHS_parallel.sh` \
This script generates gene-level iHS scores in parallel for 22 chromosomes using `iHS_results.R` (goes through each chromosome)and `iHS_scan_lite.R` (goes through each gene).

### scripts/SDS/

This directory includes scripts required to generate gene-level SDS scores. Paths are absolute with regard to the cluster, so it won't be possible to run these scripts on another machine without significant editing. The output of these scripts is included in `figures/data/SDS_results.txt.gz`. 

- `SDS_split.sh` \
This script prepares UK10k SDS scores gathered from [Field et al. 2016](https://www.science.org/doi/10.1126/science.aag0776?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed). 

- `SDS_parallel.sh` \
This script generates gene-level iHS scores in parallel for 22 chromosomes using `SDS_results.R` (goes through each chromosome)and `SDS_scan_lite.R` (goes through each gene).

### scripts/data/

This directory includes all files required to run `expression.R`, `p_cutoff_time_series_regression.R` and `random_time_series_regression.R`. 


