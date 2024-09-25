#!/bin/bash

# Define directories and input files
hicdir="path/to/hic_data"  # Set the directory path for Hi-C data
fq1="$hicdir/hic_1.fastq.gz"  # Hi-C reads R1
fq2="$hicdir/hic_2.fastq.gz"  # Hi-C reads R2
assembly="contig.asm.fasta"  # Primary contig assembly file
output_prefix="sample"  # Prefix for output files
restriction_enzyme="HindIII"  # Restriction enzyme used (HindIII or MboI)

# Step 1: Build index for the assembly
bwa index ${assembly}
samtools faidx ${assembly}

# Step 2: Align Hi-C reads using BWA
bwa aln -t 6 ${assembly} ${fq1} > ${output_prefix}_R1.sai
bwa aln -t 6 ${assembly} ${fq2} > ${output_prefix}_R2.sai
bwa sampe ${assembly} ${output_prefix}_R1.sai ${output_prefix}_R2.sai ${fq1} ${fq2} > ${output_prefix}.bwa_aln.sam

# Step 3: Filtering the SAM file
ALLHiC/scripts/PreprocessSAMs.pl ${output_prefix}.bwa_aln.sam ${assembly} ${restriction_enzyme}
ALLHiC/scripts/filterBAM_forHiC.pl ${output_prefix}.bwa_aln.REduced.paired_only.bam ${output_prefix}.clean.sam
samtools view -bt ${assembly}.fai ${output_prefix}.clean.sam > ${output_prefix}.clean.bam

# Step 4: Partition the BAM file
ALLHiC_partition -b ${output_prefix}.clean.bam -r ${assembly} -e ${restriction_enzyme} -k 21

# Step 5: Extract CLM file and count restriction sites
# For HindIII use AAGCTT; for MboI use GATC
allhic extract ${output_prefix}.clean.bam ${assembly} --RE AAGCTT

# Step 6: Optimize the assembly based on restriction enzyme counts
# This loop optimizes over several iterations (adjust $k as necessary)
k=3  # Adjust k based on specific optimization needs
for K in `seq 1 $k`;
do
    allhic optimize ${output_prefix}.clean.counts_AAGCTT.${k}g${K}.txt ${output_prefix}.clean.clm
done

# Step 7: Final assembly build
ALLHiC_build ${assembly}

# Step 8: Build index for the final assembly groups
samtools faidx groups.asm.fasta

# Create a chromosome list file for plotting
cut -f 1-2 groups.asm.fasta.fai | grep "${output_prefix}.clean.counts" > groups.chrn.list

# Step 9: Generate a Hi-C heatmap plot
ALLHiC_plot -b ${output_prefix}.clean.bam -a groups.agp -l groups.chrn.list -m 500k -o groups.hic.heatmap

