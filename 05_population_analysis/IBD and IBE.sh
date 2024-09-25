#!/bin/bash

# Define input files and directories
VCF_FILE="input_snps.vcf.gz"  # Input VCF file with SNP data
OUTPUT_DIR="IBD_IBE_analysis"  # Output directory for results
POPULATIONS_FILE="populations.txt"  # File containing population groups for FST calculations
GEOGRAPHIC_MATRIX="geographic_distance_matrix.txt"  # Matrix of geographic distances in km
ENVIRONMENTAL_MATRIX="environmental_distance_matrix.txt"  # Matrix of environmental distances

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Calculate pairwise FST values using VCFtools with a window size of 100 kb and step size of 10 kb
vcftools --gzvcf ${VCF_FILE} --weir-fst-pop ${POPULATIONS_FILE} --fst-window-size 100000 --fst-window-step 10000 --out ${OUTPUT_DIR}/fst_results

# Step 2: Convert FST values into FST/(1-FST) for Mantel test input
awk '{if ($5 != "nan") print $1, $2, $3, $4, $5 / (1 - $5)}' ${OUTPUT_DIR}/fst_results.windowed.weir.fst > ${OUTPUT_DIR}/fst_transformed.txt

# Prepare data for Mantel test (R analysis)
# For this analysis, you need three matrices: FST/(1-FST), geographic distances, and environmental distances.
# Ensure that these matrices are correctly formatted for the Mantel test in R.

echo "Preparation for IBD and IBE analysis completed. Next steps: Perform Mantel test using R and the vegan package."

# Sample R script to run Mantel test (save this as 'mantel_test.R' and run in R):
cat <<EOT > ${OUTPUT_DIR}/mantel_test.R
# Load necessary libraries
library(vegan)

# Load data matrices
fst_matrix <- as.matrix(read.table("${OUTPUT_DIR}/fst_transformed.txt", header=TRUE, row.names=1))
geo_matrix <- as.matrix(read.table("${GEOGRAPHIC_MATRIX}", header=TRUE, row.names=1))
env_matrix <- as.matrix(read.table("${ENVIRONMENTAL_MATRIX}", header=TRUE, row.names=1))

# Mantel test for IBD (FST vs. geographic distance)
mantel_result_IBD <- mantel(fst_matrix, geo_matrix, permutations = 9999, method = "pearson")
print("Mantel test for Isolation by Distance (IBD):")
print(mantel_result_IBD)

# Mantel test for IBE (FST vs. environmental distance)
mantel_result_IBE <- mantel(fst_matrix, env_matrix, permutations = 9999, method = "pearson")
print("Mantel test for Isolation by Environment (IBE):")
print(mantel_result_IBE)
EOT

echo "Run the R script for Mantel test: Rscript ${OUTPUT_DIR}/mantel_test.R"
