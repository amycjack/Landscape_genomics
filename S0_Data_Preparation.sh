#!/bin/bash 
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J snp_manipulation
#SBATCH -t 01:00:00

## Step 1: Preparing SNP Data for LFMM Analysis
### Filtering SNP Data
# We filter the SNP data by a minor allele frequency (MAF) of 0.05 and randomly select 50,000 SNPs (to reduce computational complexity) for LFMM analysis.

module load vcftools bcftools

# Load Plink and VCFtools to prepare SNP input for LFMM.
# Filter vcf by 0.05 maf.
vcftools --vcf fonio.Dexilis.vcf --maf 0.05 --recode --out fonio.Dexilis.05maf.vcf

# Randomly select 50k SNPs.
bcftools sort fonio.Dexilis.vcf -Oz -o sorted_input.vcf.gz --temp-dir tmp
bcftools index sorted_input.vcf.gz
bcftools view -H sorted_input.maf05.vcf | shuf -n 50000 | sort -k1,1 -k2,2n > random_snps.txt
bcftools view -T random_snps.txt sorted_input.maf05.vcf -Oz -o selected_snps.vcf.gz

# Convert to VCF and to genotypes envoded as 0, 1, or 2 for for homozygous, heterozygous, and other homozygous, respectively.
vcftools --gzvcf selected_snps.vcf.gz --012 --out snp
sed 's/-1/9/g' snp.012 | cut -f2- > snp.lfmm


## Step 2: Preparing SNP Data for GF Modelling

# Prepare the snp.012 that you made for lfmm and alter for GF model.
cut -f2- snp.012 | sed 's/-1/NA/g' >snp.temp
tr '\t' '_' <snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - snp.012.indv) <(echo "" | cat header - snp.temp) > snp.forR
