# Define directories and input files
hicdir = "path/to/hic_data"  # Set the directory path for Hi-C data
fq1 = f"{hicdir}/hic_1.fastq.gz"  # Hi-C reads R1
fq2 = f"{hicdir}/hic_2.fastq.gz"  # Hi-C reads R2
assembly = "contig.asm.fasta"  # Primary contig assembly file
output_prefix = "sample"  # Prefix for output files
restriction_enzyme = "HindIII"  # Restriction enzyme used (HindIII or MboI)
k = 3  # Number of iterations for optimization

# Rule all - Specifies the final outputs of the workflow
rule all:
    input:
        f"{output_prefix}.clean.bam",
        "groups.hic.heatmap"

# Step 1: Build index for the assembly
rule build_index:
    input:
        assembly=assembly
    output:
        index=f"{assembly}.bwt",
        fai=f"{assembly}.fai"
    shell:
        """
        bwa index {input.assembly}
        samtools faidx {input.assembly}
        """

# Step 2: Align Hi-C reads using BWA
rule align_hic_reads:
    input:
        assembly=assembly,
        fq1=fq1,
        fq2=fq2
    output:
        sam=f"{output_prefix}.bwa_aln.sam",
        sai1=f"{output_prefix}_R1.sai",
        sai2=f"{output_prefix}_R2.sai"
    threads: 6
    shell:
        """
        bwa aln -t {threads} {input.assembly} {input.fq1} > {output.sai1}
        bwa aln -t {threads} {input.assembly} {input.fq2} > {output.sai2}
        bwa sampe {input.assembly} {output.sai1} {output.sai2} {input.fq1} {input.fq2} > {output.sam}
        """

# Step 3: Filtering the SAM file
rule filter_sam:
    input:
        sam=f"{output_prefix}.bwa_aln.sam",
        assembly=assembly
    output:
        clean_bam=f"{output_prefix}.clean.bam"
    shell:
        """
        ALLHiC/scripts/PreprocessSAMs.pl {input.sam} {input.assembly} {restriction_enzyme}
        ALLHiC/scripts/filterBAM_forHiC.pl {output.clean_bam}.REduced.paired_only.bam {output.clean_bam}.sam
        samtools view -bt {input.assembly}.fai {output.clean_bam}.sam > {output.clean_bam}
        """

# Step 4: Partition the BAM file
rule partition_bam:
    input:
        clean_bam=f"{output_prefix}.clean.bam",
        assembly=assembly
    shell:
        """
        ALLHiC_partition -b {input.clean_bam} -r {input.assembly} -e {restriction_enzyme} -k 21
        """

# Step 5: Extract CLM file and count restriction sites
rule extract_clm:
    input:
        clean_bam=f"{output_prefix}.clean.bam",
        assembly=assembly
    output:
        clm_file=f"{output_prefix}.clean.counts_AAGCTT.clm"
    shell:
        """
        allhic extract {input.clean_bam} {input.assembly} --RE AAGCTT
        """

# Step 6: Optimize the assembly based on restriction enzyme counts
rule optimize_assembly:
    input:
        clm_file=f"{output_prefix}.clean.counts_AAGCTT.clm"
    output:
        optimized_clm=f"{output_prefix}.optimized.clm"
    shell:
        """
        for K in $(seq 1 {k});
        do
            allhic optimize {input.clm_file}.${{k}}g${{K}}.txt {input.clm_file} > {output.optimized_clm}
        done
        """

# Step 7: Final assembly build
rule build_final_assembly:
    input:
        assembly=assembly
    output:
        groups_asm="groups.asm.fasta"
    shell:
        """
        ALLHiC_build {input.assembly}
        """

# Step 8: Build index for the final assembly groups
rule build_index_final_assembly:
    input:
        groups_asm="groups.asm.fasta"
    output:
        groups_fai="groups.asm.fasta.fai"
    shell:
        """
        samtools faidx {input.groups_asm}
        """

# Step 9: Generate a Hi-C heatmap plot
rule generate_heatmap:
    input:
        clean_bam=f"{output_prefix}.clean.bam",
        fai="groups.asm.fasta.fai"
    output:
        heatmap="groups.hic.heatmap"
    shell:
        """
        cut -f 1-2 {input.fai} | grep "{output_prefix}.clean.counts" > groups.chrn.list
        ALLHiC_plot -b {input.clean_bam} -a groups.agp -l groups.chrn.list -m 500k -o {output.heatmap}
        """


