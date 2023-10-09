#!/bin/bash 

# Lin Poyraz 06/28/22
# Usage: create_dosages.sh 

module load bcftools 
module load python 

#my_dir=/project/mathilab/ppoyraz/dosage 
#my_dir=/project/mathilab/ppoyraz/1240k_dosage
my_dir=/project/mathilab/ppoyraz/dosage_updated
temp=${my_dir}/temp
brit=/project/mathilab/data/brit_data2/Imputed_hg38_TopMed
modern=/project/mathilab/data/1000g/NYGC
jti=/project/mathilab/colbranl/data/jti_snp_coordinates.txt.gz
bin=/project/mathilab/ppoyraz/time_series_predixcan/bin 
analysis=/project/mathilab/ppoyraz/time_series_predixcan/analysis
ids=/project/mathilab/colbranl/data/gtex_v8_vcfSNPs_hg38_refFix.txt.gz
mkdir -p ${my_dir}
mkdir -p ${temp}
mkdir -p ${my_dir}/dosage_files
mkdir -p ${my_dir}/1240k_dosage_files
mkdir -p ${my_dir}/missing

for chr in {1..22};do
  echo "Preparing chromosome ${chr} files"

  echo "Creating JTI SNP file..."
  #zgrep ^chr${chr}$'\t' ${jti} | cut -f1,3 | sort -nk2 | grep -v chr${chr}_ > ${temp}/jti_chr${chr}_snps

  echo "Preparing ancient british data..."
  #bsub "bcftools view -T ${temp}/jti_chr${chr}_snps -Oz ${brit}/chr${chr}.dose.vcf.gz > ${temp}/jti_brit_chr${chr}.dose.vcf.gz"

  echo "Preparing 1000g data..."
  #bsub "bcftools view -Oz -T ${temp}/jti_chr${chr}_snps -S ${my_dir}/GBR_samples ${modern}/CCDG_13607_B01_GRM_WGS_2019-02-19_chr${chr}.recalibrated_variants.vcf.gz > ${temp}/jti_1000G_GBR_chr${chr}.vcf.gz"
  
  echo "Merging 1000g and brit data..."  
  #bcftools index -f ${temp}/jti_brit_chr${chr}.dose.vcf.gz
  #bcftools index -f ${temp}/jti_1000G_GBR_chr${chr}.vcf.gz
  #bcftools merge -Ou ${temp}/jti_1000G_GBR_chr${chr}.vcf.gz ${temp}/jti_brit_chr${chr}.dose.vcf.gz | bcftools annotate -Oz -x "^INFO/AF,^FORMAT/DS,^FORMAT/GT" > ${temp}/jti_merged_chr${chr}_snps.vcf.gz

  echo "Converting to dosage format..."
  #python ${bin}/vcf2dosageprob.py ${temp}/jti_merged_chr${chr}_snps.vcf.gz > ${temp}/jti_merged_chr${chr}_dosage.txt
  #gzip ${temp}/jti_merged_chr${chr}_dosage.txt
  
  echo "Updating rsIDs..."
  # 5 1 3 4 and {ids} for 1240k db
  # 3 2 4 6 and {jti} for bestmodels db 
  #python ${bin}/update_rsIDs.py dosage 3 2 4 6 ${jti} ${my_dir}/dosage_files ${temp}/jti_merged_chr${chr}_dosage.txt.gz
  
done

Rscript ${analysis}/predixcan/prepare_predixcan_dir.R ${my_dir}/dosage_files

#gzip -f ${my_dir}/dosage_files/*.dos.txt

#Rscript ${analysis}/predixcan/prepare_predixcan_dir.R ${my_dir}/1240k_dosage_files

#gzip ${my_dir}/1240k_dosage_files/*.dos.txt

#bsub -M 24000 -e preds_%J.error -o preds_%J.out "./run_1240k_predixcan.sh predixcan_1240k_all_best ${my_dir}/1240k_dosage_files ${my_dir}"

bsub -M 24000 -e preds_%J.error -o preds_%J.out "./run_predixcan.sh Best_Models_ProteinCoding ${my_dir}/dosage_files ${my_dir}"

#bsub -M 24000 -e preds_%J.error -o preds_%J.out "./run_matched_predixcan.sh full_1240k_tissues_matched ${my_dir}/dosage_files ${my_dir}/full_matched_predictions"

#bsub -M 24000 -e preds_%J.error -o preds_%J.out "./run_matched_predixcan.sh 1240k_full_tissues_matched ${my_dir}/dosage_files ${my_dir}/1240k_matched_predictions"

#Rscript ./report_missing_snps.R ${my_dir}
