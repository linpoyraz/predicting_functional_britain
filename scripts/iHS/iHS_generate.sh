#!/bin/bash 

# Lin Poyraz 
# Adapted from Laura Colbran
# This script generates iHS scores for 1kg variants using selscan. 

module load bcftools 
module load R
module load htslib 

my_dir=/project/mathilab/ppoyraz/scan/iHS_selscan
jti_snps=/project/mathilab/ppoyraz/dosage/temp
temp=${my_dir}/temp
data=/project/mathilab/data/1000g/NYGC_phased
annots=/project/mathilab/data/1000g/NYGC
jti=/project/mathilab/colbranl/data/jti_snp_coordinates.txt.gz
bin=/project/mathilab/ppoyraz/time_series_predixcan/bin 
GBR=/project/mathilab/ppoyraz/dosage/GBR_samples
chain=/project/mathilab/ppoyraz/lift/hg38ToPanTro6.over.chain.gz

mkdir -p ${my_dir}
mkdir -p ${my_dir}/results_AA
mkdir -p ${temp}
mkdir -p ${temp}/interp_maps_AA

# ensure that ref/alt matches ancestral/derived encoding 
for chr in {1..22};do

	bsub "bcftools view ${data}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -S ${GBR} -m2 -M2 -v snps | bcftools view - -g ^miss -q 0.05 | bgzip > ${temp}/1kG_chr${chr}_GBR.vcf_not_annot.gz
	bsub "bcftools index ${temp}/1kG_chr${chr}_GBR.vcf_not_annot.gz

	bcftools annotate -Oz -a ${annots}/CCDG_13607_B01_GRM_WGS_2019-02-19_chr${chr}.recalibrated_variants.vcf.gz -c ID ${temp}/1kG_chr${chr}_GBR.vcf_not_annot.gz > ${temp}/1kG_chr${chr}_GBR.vcf.gz"

	echo "Creating bed files from vcf files"
	zcat ${temp}/1kG_chr${chr}_GBR.vcf.gz | awk -v OFS='\t' -v FS='\t' '/^[^#]/ {print $1, $2, $2, $3}' > ${temp}/chr${chr}.bed
	
	bsub "~/.local/bin/liftOver ${temp}/chr${chr}.bed ${chain} ${temp}/output_chr${chr}.bed ${temp}/unlifted_chr${chr}.bed"

	echo "Creating AA files"
	Rscript ancestral_match.R ${temp}/output_chr${chr}.bed ${temp}/unlifted_chr${chr}.bed ${temp}/1kG_chr${chr}_GBR.vcf.gz

	echo "Annotating vcf with AA information, switching ancestral/ref based on annotation"
 	# extract SNPs that have anc allele 
	awk '{print $1}' ${temp}/chr${chr}_AA.txt | grep ^rs > ${temp}/chr${chr}_AA_filter.txt

	# Change ref allele coding to ancestral and convert to vcf
	~/.local/bin/plink2 --vcf ${temp}/1kG_chr${chr}_GBR.vcf.gz --extract ${temp}/chr${chr}_AA_filter.txt --ref-allele force ${temp}/chr${chr}_AA.txt --recode vcf --out ${temp}/chr${chr}_AA_filtered 
	bsub "gzip -f ${temp}/chr${chr}_AA_filtered.vcf
	
	zcat ${temp}/chr${chr}_AA_filtered.vcf.gz | grep -v ^# | cut -f 2 > ${temp}/GBR_chr${chr}_newpos_AA.txt
   	
	./interpolate_map.R chr${chr} /project/mathilab/data/maps/1kG_hg38/GBR/chr${chr}_hg38.map ${temp}/GBR_chr${chr}_newpos_AA.txt ${temp}/interp_maps_AA/GBR_chr${chr}.map
   	
	~/.local/bin/selscan --ihs --map ${temp}/interp_maps_AA/GBR_chr${chr}.map --vcf ${temp}/chr${chr}_AA_filtered.vcf.gz --out ${my_dir}/results_AA/GBR_chr${chr}"
	
done  

~/.local/bin/norm --ihs --files ${my_dir}/results_AA/*.ihs.out

Rscript iHS_split.R 