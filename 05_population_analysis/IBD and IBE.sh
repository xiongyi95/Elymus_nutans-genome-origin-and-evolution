#!/bin/bash

mkdir -p IBD_IBE_analysis

# Step 1: Calculate pairwise FST values using VCFtools with a window size of 100 kb and step size of 10 kb
vcftools --gzvcf input_snps.vcf.gz --weir-fst-pop populations.txt --fst-window-size 100000 --fst-window-step 10000 --out IBD_IBE_analysis/esults

# Step 2: Convert FST values into FST/(1-FST) for Mantel test input
awk '{if ($5 != "nan") print $1, $2, $3, $4, $5 / (1 - $5)}' IBD_IBE_analysis/fst_results.windowed.weir.fst > IBD_IBE_analysis/fst_transformed.txt

# Prepare data for Mantel test (R analysis)
# For this analysis, you need three matrices: FST/(1-FST), geographic distances, and environmental distances.
# Ensure that these matrices are correctly formatted for the Mantel test in R.

echo "Preparation for IBD and IBE analysis completed. Next steps: Perform Mantel test using R and the vegan package."

# Sample R script to run Mantel test (save this as 'mantel_test.R' and run in R):
cat <<EOT > mantel_test.R
# Load necessary libraries
library(vegan)

# Load data matrices
fst_matrix <- as.matrix(read.table("${OUTPUT_DIR}/fst_transformed.txt", header=TRUE, row.names=1))
geo_matrix <- as.matrix(read.table("environmental_distance_matrix.txt", header=TRUE, row.names=1))
env_matrix <- as.matrix(read.table("environmental_distance_matrix.txt", header=TRUE, row.names=1))

# Mantel test for IBD (FST vs. geographic distance)
mantel_result_IBD <- mantel(fst_matrix, geo_matrix, permutations = 9999, method = "pearson")
print("Mantel test for Isolation by Distance (IBD):")
print(mantel_result_IBD)

# Mantel test for IBE (FST vs. environmental distance)
mantel_result_IBE <- mantel(fst_matrix, env_matrix, permutations = 9999, method = "pearson")
print("Mantel test for Isolation by Environment (IBE):")
print(mantel_result_IBE)
EOT

