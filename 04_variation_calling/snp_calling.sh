#!/bin/bash

# Loop through each accession ID in sample.txt
cat ${ACCESSIONS_FILE} | sed s/[[:space:]]//g | while read id; do

# Step 1: Trimming reads using Trimmomatic
java -jar ${TRIMMOMATIC_PATH}/trimmomatic-0.33.jar PE -threads 8 \
        /work/HBQB/cleandata/${id}_1.fq.gz /work/HBQB/cleandata/${id}_2.fq.gz \
        /work/HBQB/cleandata/${id}_trimmed_1.fq.gz /work/HBQB/cleandata/${id}_unpaired_1.fq.gz \
        /work/HBQB/cleandata/${id}_trimmed_2.fq.gz /work/HBQB/cleandata/${id}_unpaired_2.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Step 2: Mapping trimmed reads with BWA and removing raw files
bwa mem -K 100000000 -t 64 -M -R "@RG\tID:${id}\tPL:illumina\tLB:${id}\tSM:${id}" \
    ${REF_GENOME} /work/HBQB/cleandata/${id}_trimmed_1.fq.gz /work/HBQB/cleandata/${id}_trimmed_2.fq.gz \
    | samtools view -Sb - > /work/HBQB/bamfiles/${id}.bam
rm -f /work/HBQB/cleandata/${id}_trimmed_1.fq.gz /work/HBQB/cleandata/${id}_trimmed_2.fq.gz

# Step 3: Sorting BAM files and marking duplicates with Picard
sambamba sort --show-progress /work/HBQB/bamfiles/$id.bam -t 8 -l 3 -o /work/HBQB/bamfiles/$id.sorted.bam &&
rm -f /work/HBQB/bamfiles/$id.bam &&
sambamba markdup --show-progress -t 8 -l 3 /work/HBQB/bamfiles/$id.sorted.bam  /work/HBQB/bamfiles/$id.DEDUPED.BAM &&
sambamba index --show-progress -t 8 /work/HBQB/bamfiles/$id.DEDUPED.BAM &&
rm -f /work/HBQB/bamfiles/$id.sorted.bam /work/HBQB/bamfiles/$id.sorted.bam.bai &&

    samtools index /work/HBQB/bamfiles/${id}.deduped.bam

# Step 4: SNP calling with GATK HaplotypeCaller
gatk HaplotypeCaller -R ${REF_GENOME} -I /work/HBQB/bamfiles/${id}.deduped.bam \
-O /work/HBQB/gvcfs/${id}.g.vcf.gz -ERC GVCF

# Cleanup: Removing temporary and processed files
rm -f /work/HBQB/bamfiles/${id}.deduped.bam /work/HBQB/bamfiles/${id}.deduped.bam.bai

# Step 5: Joint Genotyping with GATK
gatk GenotypeGVCFs -R ${REF_GENOME} -V gendb://path/to/gvcfs_database \
    -O /work/HBQB/merged_population.vcf.gz
done

# Step 6: SNP filtering using criteria
gatk VariantFiltration -R ${REF_GENOME} -V /work/HBQB/merged_population.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
    --filter-name "FILTER" -O /work/HBQB/filtered_population.vcf.gz

