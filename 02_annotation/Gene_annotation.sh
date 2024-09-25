#!/bin/bash

# Create output directory
mkdir -pgene_annotation

# Step 1: Mask the genome using RepeatMasker
RepeatMasker -pa 8 -lib ${TE_LIB} -dirgene_annotation/repeatmasker_outputmasked_genome.fasta

# Step 2: RNA-seq data processing
# RNA-seq reads trimming using Trimmomatic
for fq inrna_seq_data/*.fastq.gz; do
    trimmed_fq=${OUTPUT_DIR}/trimmed_$(basename ${fq})
    trimmomatic PE -threads 8 ${fq} ${fq/.fastq.gz/_2.fastq.gz} \
        ${trimmed_fq} ${trimmed_fq/.fastq/_unpaired.fastq.gz} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Step 3: Transcriptome assembly using Trinity
Trinity --seqType fq --leftgene_annotation/trimmed_*_1.fastq.gz --rightgene_annotation/trimmed_*_2.fastq.gz \
    --CPU 16 --max_memory 100G --outputgene_annotation/trinity_output

# Step 4: Map transcripts to the genome and cluster alignments using PASA
pasa_asmbls_to_training_set.dbi --genomemasked_genome.fasta --transcriptsgene_annotation/trinity_output/Trinity.fasta \
    --CPU 16 --outputgene_annotation/pasa_output

# Step 5: Optimize gene structure using Augustus
augustus --species=e_nutans --hintsfile=./pasa_output/pasa.gff3 --extrinsicCfgFile extrinsic.ME.cfg \
    --outfile=./augustus_output.gff3masked_genome.fasta

# Step 6: Homologous annotation using GeMoMa with three related species
GeMoMa -gmasked_genome.fasta -a species1.gff -t species1.fasta -e evalue --outdirgene_annotation/gemoma_output
GeMoMa -gmasked_genome.fasta -a species2.gff -t species2.fasta -e evalue --outdirgene_annotation/gemoma_output
GeMoMa -gmasked_genome.fasta -a species3.gff -t species3.fasta -e evalue --outdirgene_annotation/gemoma_output

# Step 7: Ab initio prediction using Augustus and Helixer
augustus --species=e_nutans.masked_genome.fasta >gene_annotation/augustus_ab_initio.gff3
helixer --genomemasked_genome.fasta --outputgene_annotation/helixer_output

# Step 8: Combine all predictions into a common gene set using EVM
EVM --genomemasked_genome.fasta --transcript_alignmentsgene_annotation/pasa_output/pasa.gff3 \
    --protein_alignmentsgene_annotation/gemoma_output/*.gff3 --augustus_predictionsgene_annotation/augustus_output.gff3 \
    --helixer_predictionsgene_annotation/helixer_output.gff3 --weightsevm_weights.txt \
    --output_filegene_annotation/evm_combined.gff3

# Step 9: Annotation of non-coding RNAs (tRNA, rRNA, snRNA, miRNA)
# Predict tRNA using tRNAscan-SE
tRNAscan-SE -ogene_annotation/tRNAscan_output.gff3 -fgene_annotation/tRNAscan_stats.txtmasked_genome.fasta

# Predict rRNA using Barrnap
barrnap --kingdom euk --quietmasked_genome.fasta >gene_annotation/barrnap_output.gff3

# Predict snRNA and miRNA using Infernal and Rfam
cmsearch --cpu 16 --tbloutgene_annotation/cmsearch_output.tblout Rfam.cmmasked_genome.fasta

echo "Gene and non-coding RNA annotation completed. Results are saved ingene_annotation."
