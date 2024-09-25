#!/bin/bash

# Step 1: Principal Component Analysis (PCA) using PLINK
plink --vcf snps.vcf.gz --make-bed --out analysis_results/pca_input
plink --bfile analysis_results/pca_input --pca --out analysis_results/pca_results

# Step 2: Construct a neighbor-joining tree using FastTree
# Convert VCF to FASTA for tree construction
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' snps.vcf.gz | awk '{print ">"$1"_"$2"\n"$3$4}' > analysis_results/snps.fasta
fasttree -nt analysis_results/snps.fasta > analysis_results/nj_tree.tree

# Step 3: Annotate SNPs using SnpEff
java -Xmx4g -jar snpEff.jar -v reference_genome snps.vcf.gz > analysis_results/annotated_snps.vcf

# Step 4: Analyze linkage disequilibrium (LD) with PopLDdecay
# Convert VCF to PopLDdecay input format and analyze LD
VCF2Dis -InVCF snps.vcf.gz -OutPut analysis_results/poplddecay_input
PopLDdecay -InVCF snps.vcf.gz -OutStat analysis_results/poplddecay_stats -MaxDist 100 -MinMAF 0.05 -Miss 0.1 -WinSize 100000 -StepSize 10000

# Step 5: Calculate nucleotide diversity (π) using VCFtools
# Split the VCF file by population and calculate π for each population
for POP in $(cat populations_list.txt); do
    vcftools --gzvcf snps.vcf.gz --keep ${POP} --window-pi 100000 --out analysis_results/nucleotide_diversity_${POP}
done

# Step 6: Calculate mean π values across the whole genome for each population
for POP in $(cat populations_list.txt); do
    awk '{sum += $4; n++} END {if (n > 0) print "Mean π for " POP ":", sum / n; else print "No data available for " POP}' analysis_results/nucleotide_diversity_${POP}.windowed.pi
done


