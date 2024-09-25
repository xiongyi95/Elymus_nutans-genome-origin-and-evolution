#!/bin/bash
# Step 1: Identify single-copy orthologs using OrthoFinder
orthofinder -f $GENOME_DIR -t 96 -o $OUTPUT_DIR/orthofinder_output

# Step 2: Extract coding sequences (CDS) of single-copy orthologs
cd $OUTPUT_DIR/orthofinder_output/Results*/Orthogroups
for og in $(cat SingleCopyOrthogroups.txt); do
    grep -f <(echo "${og}") Orthogroups.txt | cut -f 2- > $OUTPUT_DIR/cds_extraction/${og}.cds
done

# Step 3: Merge aligned protein and CDS alignments to construct a supermatrix
for cds in $OUTPUT_DIR/cds_extraction/*.cds; do
    iqtree2 -s ${cds} -m MFP -bb 1000 -nt 16 -o $OUTPUT_DIR/phylogenetic_tree/$(basename ${cds} .cds).treefile
done

# Step 4: Extract 4DTv sites and estimate transversion rates
# extract_4dtv_sites.py --input supermatrix.aligned.fasta --output $OUTPUT_DIR/4dtv_sites.txt

# Step 5: Construct evolutionary tree using IQ-TREE and ASTRAL
java -jar $ASTRAL_PATH -i $OUTPUT_DIR/phylogenetic_tree/*.treefile -o $OUTPUT_DIR/phylogenetic_tree/supermatrix_astral.tree

# Step 6: Estimate divergence times using MCMCTree from PAML
cd $OUTPUT_DIR/divergence_time_estimation
mcmctree mcmctree.ctl

# Step 7: Analyze the evolutionary position of E. nutans using WGDI with the -at parameter
WGDI -g $GENOME_DIR -a $CDS_DIR -t subgenomes_list.txt -o $OUTPUT_DIR/wgdi_output

# Step 8: Align the high-quality collinear genes using Muscle
muscle -align $OUTPUT_DIR/wgdi_output/collinear_genes.fasta -output $OUTPUT_DIR/aligned_genes.afa

# Step 9: Construct a phylogenetic tree with IQ-TREE
iqtree2 -s $OUTPUT_DIR/aligned_genes.afa -m MFP -bb 1000 -nt 16 -o $OUTPUT_DIR/phylogenetic_tree/elymus_tree.treefile

# Step 10: Merge multiple phylogenetic trees using ASTRAL
java -jar $ASTRAL_PATH -i $OUTPUT_DIR/phylogenetic_tree/*.treefile -o $OUTPUT_DIR/phylogenetic_tree/merged_astral.tree

# Step 11: Construct the chloroplast phylogenetic tree
for reads in $RAW_READS_DIR/*; do
    get_organelle_from_reads.py -1 ${reads}/reads_R1.fastq.gz -2 ${reads}/reads_R2.fastq.gz -t 16 -o $OUTPUT_DIR/chloroplast_assembly
done

# Combine assembled chloroplast genomes with downloaded NCBI chloroplast genomes
cat $OUTPUT_DIR/chloroplast_assembly/*.fasta $NCBI_CPGENOMES/*.fasta > $OUTPUT_DIR/chloroplast_tree/all_cp_genomes.fasta

# Construct the chloroplast phylogenetic tree using IQ-TREE
iqtree2 -s $OUTPUT_DIR/chloroplast_tree/all_cp_genomes.fasta -m MFP -bb 1000 -nt 16 -o $OUTPUT_DIR/chloroplast_tree/cp_phylogenetic_tree.treefile
