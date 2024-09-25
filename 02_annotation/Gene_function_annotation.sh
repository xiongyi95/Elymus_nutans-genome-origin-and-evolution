#!/bin/bash

THREADS=96  # Number of threads to use for computations

# Create output directory
mkdir -p .

# Step 1: Search against SwissProt database
blastp -query protein_coding_genes.fasta -db swissprot -out ./swissprot_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 2: Search against non-redundant protein (NR) database
blastp -query protein_coding_genes.fasta -db nr -out ./nr_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 3: KEGG annotation using KAAS
# This requires submitting to the KAAS web service, so saving sequences to a file for submission
cp protein_coding_genes.fasta ./kegg_input.fasta
# For automated submission, consider setting up a script or contacting the service with appropriate credentials

# Step 4: KOG annotation
blastp -query protein_coding_genes.fasta -db kog -out ./kog_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 5: Search against TrEMBL database
blastp -query protein_coding_genes.fasta -db trembl -out ./trembl_annotation.txt -evalue 1e-5 -num_threads ${THREADS} -outfmt 6

# Step 6: InterProScan to identify motifs, domains, and GO annotations
interproscan.sh -i protein_coding_genes.fasta -f tsv -dp -T ./interproscan_temp -o ./interproscan_annotation.tsv -cpu ${THREADS}

# Step 7: Extract GO annotations from InterProScan results
awk -F'\t' '{if ($12 != "-") print $1"\t"$12}' ./interproscan_annotation.tsv > ./go_annotations_from_interproscan.tsv

# Step 8: Predict protein functions using DeepGO
deepgoplus --data-root <path_to_data_folder> --in-file protein_coding_genes.fasta 

