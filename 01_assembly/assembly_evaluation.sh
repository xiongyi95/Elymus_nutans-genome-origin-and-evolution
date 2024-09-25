# Set input and output paths
ASSEMBLY = "path/to/assembly.fasta"  # Path to the HiFi assembly result file
HIFI_READS = "path/to/hifi_reads.fastq"  # Path to HiFi sequencing data file
HIC_READS_1 = "path/to/hic_reads_R1.fastq"  # Path to Hi-C sequencing data file (R1)
HIC_READS_2 = "path/to/hic_reads_R2.fastq"  # Path to Hi-C sequencing data file (R2)
BUSCO_DB = "path/to/busco_db"  # Path to BUSCO database
OUTPUT_DIR = "assembly_evaluation"  # Output directory
CHROMOSOME_PAIR_NUM = 14  # Set according to the chromosome number of E. nutans

# Rule all - Specifies the final outputs of the workflow
rule all:
    input:
        f"{OUTPUT_DIR}/continuity_evaluation.txt",
        f"{OUTPUT_DIR}/accuracy_evaluation.txt",
        f"{OUTPUT_DIR}/busco_results/short_summary.txt"

# Step 1: Continuity evaluation - Calculate Contig N50 and CC Ratio
rule calc_continuity:
    input:
        assembly=ASSEMBLY
    output:
        stats=f"{OUTPUT_DIR}/contig_stats.txt",
        evaluation=f"{OUTPUT_DIR}/continuity_evaluation.txt"
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        seqkit stats {input.assembly} > {output.stats}
        CONTIG_COUNT=$(grep -c ">" {input.assembly})
        CC_RATIO=$(echo "$CONTIG_COUNT / {CHROMOSOME_PAIR_NUM}" | bc -l)
        echo "CC Ratio: $CC_RATIO" >> {output.evaluation}
        """

# Step 2: Accuracy evaluation - Map HiFi reads to the assembly
rule map_hifi_reads:
    input:
        assembly=ASSEMBLY,
        hifi_reads=HIFI_READS
    output:
        bam=f"{OUTPUT_DIR}/mapped_hifi_reads.bam",
        bai=f"{OUTPUT_DIR}/mapped_hifi_reads.bam.bai"
    shell:
        """
        bwa index {input.assembly}
        bwa mem {input.assembly} {input.hifi_reads} | samtools view -b - | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule evaluate_accuracy:
    input:
        bam=f"{OUTPUT_DIR}/mapped_hifi_reads.bam"
    output:
        eval=f"{OUTPUT_DIR}/accuracy_evaluation.txt",
        variants=f"{OUTPUT_DIR}/raw_variants.vcf"
    shell:
        """
        MAPPING_RATE=$(samtools flagstat {input.bam} | grep "mapped (" | awk '{{print $5}}')
        GENOMIC_COVERAGE=$(samtools depth {input.bam} | awk '{{sum+=$3}} END {{print sum/NR}}')
        MEAN_DEPTH=$(samtools depth {input.bam} | awk '{{sum+=$3}} END {{print sum/NR}}')
        echo -e "Mapping Rate: $MAPPING_RATE\nGenomic Coverage: $GENOMIC_COVERAGE\nMean Depth: $MEAN_DEPTH" > {output.eval}
        gatk --java-options "-Xmx4g" HaplotypeCaller -R {ASSEMBLY} -I {input.bam} -O {output.variants}
        SNP_COUNT=$(grep -v "^#" {output.variants} | wc -l)
        echo "SNP Count: $SNP_COUNT" >> {output.eval}
        """

# Step 3: Completeness evaluation - Use BUSCO
rule busco_completeness:
    input:
        assembly=ASSEMBLY
    output:
        busco_dir=directory(f"{OUTPUT_DIR}/busco_results")
    shell:
        """
        busco -i {input.assembly} -l {BUSCO_DB} -o {output.busco_dir} -m genome
        """
