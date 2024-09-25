#!/bin/bash

# Set paths for input and output directories
HIFI_READS="path/to/hifi_reads.fastq"  # Path to Pacbio HiFi reads
OUTPUT_DIR="assembly_output"  # Output directory for assembly results
ASSEMBLY_PREFIX="assembly"  # Prefix for assembly output files

# Create output directory
mkdir -p ${OUTPUT_DIR}

#  De novo assembly with Hifiasm
hifiasm -o ${OUTPUT_DIR}/${ASSEMBLY_PREFIX} -t 16 ${HIFI_READS}

# Convert GFA to FASTA (primary contigs)
awk '$1 ~/S/{print ">"$2"\n"$3}' ${OUTPUT_DIR}/${ASSEMBLY_PREFIX}.bp.hap1.p_ctg.gfa > ${OUTPUT_DIR}/${ASSEMBLY_PREFIX}.primary_contigs.fasta
