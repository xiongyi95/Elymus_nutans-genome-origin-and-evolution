#!/bin/bash

# Define input paths and output directories
E_NUTANS_PROTEINS="path/to/elymus_nutans_proteins.fasta"  # Path to the protein sequences of E. nutans
RICE_PROTEINS="path/to/rice_proteins.fasta"  # Path to the protein sequences of rice
OUTPUT_DIR="karyotype_analysis"  # Output directory for analysis results
WGDI_PATH="/path/to/wgdi"  # Path to the WGDI software
PAML_PATH="/path/to/paml"  # Path to PAML package (YN00 for Ka/Ks calculation)

# Step 1: Merge protein sequences of rice and E. nutans and perform alignment using diamond_blastp.pl
cat ${RICE_PROTEINS} ${E_NUTANS_PROTEINS} > ${OUTPUT_DIR}/merged_proteins.fasta
${WGDI_PATH}/scripts/diamond_blastp.pl -i ${OUTPUT_DIR}/merged_proteins.fasta -o ${OUTPUT_DIR}/blastp_results.txt

# Step 2: Extract collinear genes using the -icl parameter
WGDI -i ${OUTPUT_DIR}/blastp_results.txt -icl -o ${OUTPUT_DIR}/collinear_genes.txt

# Step 3: Estimate Ka and Ks substitution rates using the Nei-Gojobori method in YN00 (PAML)
# This step assumes the collinear gene pairs are prepared for YN00 input
# Extract sequences and prepare input for YN00 (requires sequence alignment preprocessing)
WGDI -i ${OUTPUT_DIR}/collinear_genes.txt -ks -o ${OUTPUT_DIR}/ka_ks_results.txt
cd ${OUTPUT_DIR}
${PAML_PATH}/bin/yn00 yn00.ctl  # Ensure yn00.ctl is properly configured with correct input files

# Step 4: Integrate results from steps 2 and 3 using the -bi parameter
WGDI -i ${OUTPUT_DIR}/collinear_genes.txt -bi ${OUTPUT_DIR}/ka_ks_results.txt -o ${OUTPUT_DIR}/integrated_results.txt

# Step 5: Perform polyploid classification using the -pc parameter
WGDI -i ${OUTPUT_DIR}/integrated_results.txt -pc -o ${OUTPUT_DIR}/polyploid_classification.txt

# Step 6: Map the rice karyotype to the chromosomes of E. nutans using the -km and -m parameters
WGDI -i ${OUTPUT_DIR}/integrated_results.txt -km -m -o ${OUTPUT_DIR}/karyotype_mapping_results.txt

