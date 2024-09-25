#!/bin/bash

# Set input paths and parameters
GENOME="path/to/masked_genome.fasta"  # Masked genome assembly file from RepeatMasker
TE_LIB="path/to/te_library.fasta"  # Transposable element library created previously
RNA_SEQ_DIR="path/to/rna_seq_data"  # Directory containing RNA-seq data
OUTPUT_DIR="gene_annotation"  # Output directory for the results
SPECIES="E_nutans"  # Species name for Augustus
AUGUSTUS_SPECIES="e_nutans"  # Species model for Augustus
EVM_WEIGHTS="path/to/evm_weights.txt"  # Weights file for EVM

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Mask the genome using RepeatMasker
RepeatMasker -pa 8 -lib ${TE_LIB} -dir ${OUTPUT_DIR}/repeatmasker_output ${GENOME}

# Step 2: RNA-seq data processing
# RNA-seq reads trimming using Trimmomatic
for fq in ${RNA_SEQ_DIR}/*.fastq.gz; do
    trimmed_fq=${OUTPUT_DIR}/trimmed_$(basename ${fq})
    trimmomatic PE -threads 8 ${fq} ${fq/.fastq.gz/_2.fastq.gz} \
        ${trimmed_fq} ${trimmed_fq/.fastq/_unpaired.fastq.gz} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Step 3: Transcriptome assembly using Trinity
Trinity --seqType fq --left ${OUTPUT_DIR}/trimmed_*_1.fastq.gz --right ${OUTPUT_DIR}/trimmed_*_2.fastq.gz \
    --CPU 16 --max_memory 100G --output ${OUTPUT_DIR}/trinity_output

# Step 4: Map transcripts to the genome and cluster alignments using PASA
pasa_asmbls_to_training_set.dbi --genome ${GENOME} --transcripts ${OUTPUT_DIR}/trinity_output/Trinity.fasta \
    --CPU 16 --output ${OUTPUT_DIR}/pasa_output

# Step 5: Optimize gene structure using Augustus
augustus --species=${AUGUSTUS_SPECIES} --hintsfile=${OUTPUT_DIR}/pasa_output/pasa.gff3 --extrinsicCfgFile extrinsic.ME.cfg \
    --outfile=${OUTPUT_DIR}/augustus_output.gff3 ${GENOME}

# Step 6: Homologous annotation using GeMoMa with three related species
GeMoMa -g ${GENOME} -a species1.gff -t species1.fasta -e evalue --outdir ${OUTPUT_DIR}/gemoma_output
GeMoMa -g ${GENOME} -a species2.gff -t species2.fasta -e evalue --outdir ${OUTPUT_DIR}/gemoma_output
GeMoMa -g ${GENOME} -a species3.gff -t species3.fasta -e evalue --outdir ${OUTPUT_DIR}/gemoma_output

# Step 7: Ab initio prediction using Augustus and Helixer
augustus --species=${AUGUSTUS_SPECIES} ${GENOME} > ${OUTPUT_DIR}/augustus_ab_initio.gff3
helixer --genome ${GENOME} --output ${OUTPUT_DIR}/helixer_output

# Step 8: Combine all predictions into a common gene set using EVM
EVM --genome ${GENOME} --transcript_alignments ${OUTPUT_DIR}/pasa_output/pasa.gff3 \
    --protein_alignments ${OUTPUT_DIR}/gemoma_output/*.gff3 --augustus_predictions ${OUTPUT_DIR}/augustus_output.gff3 \
    --helixer_predictions ${OUTPUT_DIR}/helixer_output.gff3 --weights ${EVM_WEIGHTS} \
    --output_file ${OUTPUT_DIR}/evm_combined.gff3

# Step 9: Annotation of non-coding RNAs (tRNA, rRNA, snRNA, miRNA)
# Predict tRNA using tRNAscan-SE
tRNAscan-SE -o ${OUTPUT_DIR}/tRNAscan_output.gff3 -f ${OUTPUT_DIR}/tRNAscan_stats.txt ${GENOME}

# Predict rRNA using Barrnap
barrnap --kingdom euk --quiet ${GENOME} > ${OUTPUT_DIR}/barrnap_output.gff3

# Predict snRNA and miRNA using Infernal and Rfam
cmsearch --cpu 16 --tblout ${OUTPUT_DIR}/cmsearch_output.tblout Rfam.cm ${GENOME}

echo "Gene and non-coding RNA annotation completed. Results are saved in ${OUTPUT_DIR}."
