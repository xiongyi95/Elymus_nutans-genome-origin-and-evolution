#!/bin/bash

# Define paths for input and output directories
PHENO_FILE="path/to/phenotype_data.txt"  # Path to phenotype data file with BLUP values
VCF_FILE="path/to/snp_data.vcf.gz"  # Path to high-quality SNP data
OUTPUT_DIR="gwas_results"  # Output directory for GWAS results
POP_STRUCTURE_DIR="pop_structure"  # Output directory for population structure analysis
CM_PLOT_SCRIPT="path/to/cmplot_script.R"  # Path to R script for generating Manhattan plots

# Create output directories
mkdir -p ${OUTPUT_DIR} ${POP_STRUCTURE_DIR}

# Step 1: Population structure analysis using PLINK for Q matrix and GEMMA for K matrix
# Convert VCF to PLINK format for PCA
plink --vcf ${VCF_FILE} --make-bed --out ${POP_STRUCTURE_DIR}/pop_structure
# Calculate PCA to infer population structure (Q matrix)
plink --bfile ${POP_STRUCTURE_DIR}/pop_structure --pca --out ${POP_STRUCTURE_DIR}/pca_results
# Calculate the kinship matrix (K matrix) using GEMMA
gemma -bfile ${POP_STRUCTURE_DIR}/pop_structure -gk 2 -o kinship_matrix

# Step 2: BLUP Calculation
# Example of BLUP calculation using R (run this code separately if needed):
# Rscript -e "
library(lme4);
pheno_data <- read.table('${PHENO_FILE}', header = TRUE);
blup_model <- lmer(trait ~ (1|accession) + (1|year), data = pheno_data);
blup_values <- ranef(blup_model)$accession;
write.table(blup_values, '${OUTPUT_DIR}/blup_values.txt', row.names = TRUE, col.names = FALSE, sep = '\t');
"

# Step 3: GWAS Analysis using Mixed Linear Model (MLM)
# Prepare input for GEMMA with Q and K matrices
gemma -bfile ${POP_STRUCTURE_DIR}/pop_structure -p ${PHENO_FILE} -k output/kinship_matrix.cXX.txt -c ${POP_STRUCTURE_DIR}/pca_results.eigenvec -lmm 4 -o gwas_results

# Step 4: Filter significant SNPs with the cutoff P-value (0.1/number of SNPs)
NUM_SNPS=$(zcat ${VCF_FILE} | grep -v "^#" | wc -l)
P_CUTOFF=$(echo "0.1 / ${NUM_SNPS}" | bc -l)
awk -v pval=${P_CUTOFF} '$12 < pval {print $0}' ${OUTPUT_DIR}/gwas_results.assoc.txt > ${OUTPUT_DIR}/significant_snps.txt

# Step 5: Generate Manhattan plot using CMplot R package
# Rscript -e "
library(CMplot);
data <- read.table('${OUTPUT_DIR}/gwas_results.assoc.txt', header = TRUE);
CMplot(data, plot.type = 'm', LOG10 = TRUE, threshold = c(${P_CUTOFF}), threshold.col = 'red', file.output = TRUE, main = 'Manhattan Plot');
"


