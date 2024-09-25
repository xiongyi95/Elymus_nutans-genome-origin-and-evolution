#!/bin/bash

# Define input and output paths
GENES="path/to/protein_coding_genes.fasta"  # Path to protein-coding gene sequences in FASTA format
OUTPUT_DIR="function_annotation"  # Output directory for annotation results
THREADS=96  # Number of threads to use for computations

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Search against SwissProt database
blastp -query ${GENES} -db swissprot -out ${OUTPUT_DIR}/swissprot_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 2: Search against non-redundant protein (NR) database
blastp -query ${GENES} -db nr -out ${OUTPUT_DIR}/nr_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 3: KEGG annotation using KAAS
# This requires submitting to the KAAS web service, so saving sequences to a file for submission
cp ${GENES} ${OUTPUT_DIR}/kegg_input.fasta
# For automated submission, consider setting up a script or contacting the service with appropriate credentials

# Step 4: KOG annotation
blastp -query ${GENES} -db kog -out ${OUTPUT_DIR}/kog_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 5: Search against TrEMBL database
blastp -query ${GENES} -db trembl -out ${OUTPUT_DIR}/trembl_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 6: InterProScan to identify motifs, domains, and GO annotations
interproscan.sh -i ${GENES} -f tsv -dp -T ${OUTPUT_DIR}/interproscan_temp -o ${OUTPUT_DIR}/interproscan_annotation.tsv -cpu ${THREADS}

# Step 7: Extract GO annotations from InterProScan results
awk -F'\t' '{if ($12 != "-") print $1"\t"$12}' ${OUTPUT_DIR}/interproscan_annotation.tsv > ${OUTPUT_DIR}/go_annotations_from_interproscan.tsv

# Step 8: Predict protein functions using DeepGO
deepgoplus --data-root <path_to_data_folder> --in-file ${GENES} 

