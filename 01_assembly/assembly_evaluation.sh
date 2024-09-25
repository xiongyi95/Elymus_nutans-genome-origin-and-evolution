#!/bin/bash

# Set input and output paths
ASSEMBLY="path/to/assembly.fasta"  # Path to the HiFi assembly result file
HIFI_READS="path/to/hifi_reads.fastq"  # Path to HiFi sequencing data file
HIC_READS_1="path/to/hic_reads_R1.fastq"  # Path to Hi-C sequencing data file (R1)
HIC_READS_2="path/to/hic_reads_R2.fastq"  # Path to Hi-C sequencing data file (R2)
BUSCO_DB="path/to/busco_db"  # Path to BUSCO database
OUTPUT_DIR="assembly_evaluation"  # Output directory

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Continuity evaluation - Calculate Contig N50 and CC Ratio
# Use seqkit to calculate N50
seqkit stats ${ASSEMBLY} > ${OUTPUT_DIR}/contig_stats.txt
# Calculate CC Ratio
CONTIG_COUNT=$(grep -c ">" ${ASSEMBLY})
CHROMOSOME_PAIR_NUM=14  # Set according to the chromosome number of E. nutans
CC_RATIO=$(echo "${CONTIG_COUNT} / ${CHROMOSOME_PAIR_NUM}" | bc -l)
echo "CC Ratio: ${CC_RATIO}" >> ${OUTPUT_DIR}/continuity_evaluation.txt

# Step 2: Accuracy evaluation - Map HiFi reads to the assembly (used for assessment only; HiFi reads were used for assembly)
bwa index ${ASSEMBLY}
bwa mem ${ASSEMBLY} ${HIFI_READS} | samtools view -b - | samtools sort -o ${OUTPUT_DIR}/mapped_hifi_reads.bam
samtools index ${OUTPUT_DIR}/mapped_hifi_reads.bam
# Calculate mapping rate, genomic coverage, and mean depth
MAPPING_RATE=$(samtools flagstat ${OUTPUT_DIR}/mapped_hifi_reads.bam | grep "mapped (" | awk '{print $5}')
GENOMIC_COVERAGE=$(samtools depth ${OUTPUT_DIR}/mapped_hifi_reads.bam | awk '{sum+=$3} END {print sum/NR}')
MEAN_DEPTH=$(samtools depth ${OUTPUT_DIR}/mapped_hifi_reads.bam | awk '{sum+=$3} END {print sum/NR}')
echo -e "Mapping Rate: ${MAPPING_RATE}\nGenomic Coverage: ${GENOMIC_COVERAGE}\nMean Depth: ${MEAN_DEPTH}" >> ${OUTPUT_DIR}/accuracy_evaluation.txt

# Use GATK HaplotypeCaller to evaluate base-level accuracy
gatk --java-options "-Xmx4g" HaplotypeCaller -R ${ASSEMBLY} -I ${OUTPUT_DIR}/mapped_hifi_reads.bam -O ${OUTPUT_DIR}/raw_variants.vcf
# Count SNPs to evaluate the genome accuracy
SNP_COUNT=$(grep -v "^#" ${OUTPUT_DIR}/raw_variants.vcf | wc -l)
echo "SNP Count: ${SNP_COUNT}" >> ${OUTPUT_DIR}/accuracy_evaluation.txt

# Step 3: Completeness evaluation - Use BUSCO to assess the completeness of the gene space
busco -i ${ASSEMBLY} -l ${BUSCO_DB} -o ${OUTPUT_DIR}/busco_results -m genome

# Note: Hi-C data is primarily used for scaffolding the assembly, not directly for accuracy and completeness evaluation.

echo "Assembly quality evaluation completed. Results are saved in the ${OUTPUT_DIR} directory."
