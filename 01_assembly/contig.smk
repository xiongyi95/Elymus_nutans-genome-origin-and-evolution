# Set paths for input and output directories
HIFI_READS = "path/to/hifi_reads.fastq"  # Path to Pacbio HiFi reads
OUTPUT_DIR = "assembly_output"  # Output directory for assembly results
ASSEMBLY_PREFIX = "assembly"  # Prefix for assembly output files

# Rule all - Specifies the final outputs of the workflow
rule all:
    input:
        f"{OUTPUT_DIR}/{ASSEMBLY_PREFIX}.primary_contigs.fasta"

# Step 1: De novo assembly with Hifiasm
rule hifiasm_assembly:
    input:
        hifi_reads=HIFI_READS
    output:
        gfa=f"{OUTPUT_DIR}/{ASSEMBLY_PREFIX}.bp.hap1.p_ctg.gfa"
    params:
        output_prefix=f"{OUTPUT_DIR}/{ASSEMBLY_PREFIX}"
    threads: 16
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        hifiasm -o {params.output_prefix} -t {threads} {input.hifi_reads}
        """

# Step 2: Convert GFA to FASTA (primary contigs)
rule gfa_to_fasta:
    input:
        gfa=f"{OUTPUT_DIR}/{ASSEMBLY_PREFIX}.bp.hap1.p_ctg.gfa"
    output:
        fasta=f"{OUTPUT_DIR}/{ASSEMBLY_PREFIX}.primary_contigs.fasta"
    shell:
        """
        awk '$1 ~/S/{{print ">"$2"\\n"$3}}' {input.gfa} > {output.fasta}
        """
