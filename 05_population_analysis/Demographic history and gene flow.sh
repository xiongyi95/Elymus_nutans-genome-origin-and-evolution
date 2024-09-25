#!/bin/bash

# Set paths for input and output directories
RESEQ_DATA="path/to/resequencing_data"  # Path to resequencing data
SNP_DATA="path/to/snp_data.vcf"  # Path to SNP data file
OUTGROUP="path/to/outgroup_genome.fasta"  # Outgroup genome for Dsuite analysis
OUTPUT_DIR="demographic_analysis"  # Output directory
MUTATION_RATE=3.75e-8  # Mutation rate per site per year
GENERATION_TIME=1  # Generation time in years

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: PSMC analysis to estimate historical changes in Ne
for SAMPLE in $(ls ${RESEQ_DATA}/*.bam); do
    SAMPLE_NAME=$(basename ${SAMPLE} .bam)
    samtools mpileup -uf ${OUTGROUP} ${SAMPLE} | bcftools call -c - | vcfutils.pl vcf2fq > ${OUTPUT_DIR}/${SAMPLE_NAME}.fq
    fq2psmcfa -q20 ${OUTPUT_DIR}/${SAMPLE_NAME}.fq > ${OUTPUT_DIR}/${SAMPLE_NAME}.psmcfa
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${OUTPUT_DIR}/${SAMPLE_NAME}.psmc ${OUTPUT_DIR}/${SAMPLE_NAME}.psmcfa

    # Bootstrap estimates
    for i in $(seq 100); do
        splitfa ${OUTPUT_DIR}/${SAMPLE_NAME}.psmcfa > ${OUTPUT_DIR}/${SAMPLE_NAME}_split.psmcfa
        psmc -b -N25 -t15 -r5 -p "4+25*2+4+6" -o ${OUTPUT_DIR}/bootstrap/${SAMPLE_NAME}_bootstrap${i}.psmc ${OUTPUT_DIR}/${SAMPLE_NAME}_split.psmcfa
    done
done

# Plot PSMC results
for PSMC_RESULT in ${OUTPUT_DIR}/*.psmc; do
    psmc_plot.pl -u ${MUTATION_RATE} -g ${GENERATION_TIME} ${OUTPUT_DIR}/$(basename ${PSMC_RESULT} .psmc) ${PSMC_RESULT}
done

# Step 2: GADMA analysis using dadi engine
python -m gadma --input_vcf ${SNP_DATA} --outdir ${OUTPUT_DIR}/gadma_output --engine dadi --n_clusters 3

# Step 3: Model selection using DIYABC-RF
mkdir -p ${OUTPUT_DIR}/diyabc_rf
for MODEL in {1..16}; do
    diyabc -p ${OUTPUT_DIR}/diyabc_rf/model${MODEL} -o ${OUTPUT_DIR}/diyabc_rf/model${MODEL}_results -n 10000 -b 0 -r
done

# Apply Random Forest algorithm to select the best model
Rscript -e "library('DIYABC'); runRFanalysis(path='${OUTPUT_DIR}/diyabc_rf', models=16, output='${OUTPUT_DIR}/best_model.txt')"

# Visualize the optimal model using DemesDraw
python -m demesdraw.utils.plot_model --input ${OUTPUT_DIR}/best_model.txt --output ${OUTPUT_DIR}/optimal_model_plot.png

# Step 4: Investigate gene flow using Dsuite to compute D-statistics
Dsuite Dtrios ${SNP_DATA} ${OUTGROUP} -o ${OUTPUT_DIR}/Dsuite_results
Dsuite Dinvestigate ${SNP_DATA} ${OUTGROUP} -o ${OUTPUT_DIR}/D_f_statistics

# Step 5: Estimate historic gene flow using Migrate-n
migrate-n -i ${OUTPUT_DIR}/migrate_input/migrate_input.nex -o ${OUTPUT_DIR}/migrate_output/migrate_results -b 100000 -m -t 100000

echo "Demographic history and gene flow analysis completed. Results are saved in ${OUTPUT_DIR}."
