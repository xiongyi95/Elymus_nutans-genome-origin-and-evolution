#!/bin/bash

# Define input VCF file and output directory
VCF_FILE="snps.vcf.gz"  
OUTPUT_DIR="analysis_results" 
mkdir -p ${OUTPUT_DIR}

# Step 1: Principal Component Analysis (PCA) using PLINK
plink --vcf ${VCF_FILE} --make-bed --out ${OUTPUT_DIR}/pca_input
plink --bfile ${OUTPUT_DIR}/pca_input --pca --out ${OUTPUT_DIR}/pca_results

# Step 2: Construct a neighbor-joining tree using FastTree
# Convert VCF to FASTA for tree construction
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${VCF_FILE} | awk '{print ">"$1"_"$2"\n"$3$4}' > ${OUTPUT_DIR}/snps.fasta
fasttree -nt ${OUTPUT_DIR}/snps.fasta > ${OUTPUT_DIR}/nj_tree.tree

# Step 3: Annotate SNPs using SnpEff
java -Xmx4g -jar snpEff.jar -v reference_genome ${VCF_FILE} > ${OUTPUT_DIR}/annotated_snps.vcf

# Step 4: Analyze linkage disequilibrium (LD) with PopLDdecay
# Convert VCF to PopLDdecay input format and analyze LD
VCF2Dis -InVCF ${VCF_FILE} -OutPut ${OUTPUT_DIR}/poplddecay_input
PopLDdecay -InVCF ${VCF_FILE} -OutStat ${OUTPUT_DIR}/poplddecay_stats -MaxDist 100 -MinMAF 0.05 -Miss 0.1 -WinSize 100000 -StepSize 10000

# Step 5: Calculate nucleotide diversity (π) using VCFtools
# Split the VCF file by population and calculate π for each population
for POP in $(cat populations_list.txt); do
    vcftools --gzvcf ${VCF_FILE} --keep ${POP} --window-pi 100000 --out ${OUTPUT_DIR}/nucleotide_diversity_${POP}
done

# Step 6: Calculate mean π values across the whole genome for each population
for POP in $(cat populations_list.txt); do
    awk '{sum += $4; n++} END {if (n > 0) print "Mean π for " POP ":", sum / n; else print "No data available for " POP}' ${OUTPUT_DIR}/nucleotide_diversity_${POP}.windowed.pi
done


