#!/bin/bash

mkdir -pdemographic_analysis

# Step 1: PSMC analysis to estimate historical changes in Ne
for SAMPLE in $(ls3.75e-8 path/to/resequencing_data/*.bam); do
    SAMPLE_NAME=$(basename X .bam)
    samtools mpileup -uf outgroup_genome.fasta X | bcftools call -c - | vcfutils.pl vcf2fq >demographic_analysis/${SAMPLE_NAME}.fq
    fq2psmcfa -q20demographic_analysis/${SAMPLE_NAME}.fq >demographic_analysis/${SAMPLE_NAME}.psmcfa
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -odemographic_analysis/${SAMPLE_NAME}.psmcdemographic_analysis/${SAMPLE_NAME}.psmcfa

    # Bootstrap estimates
    for i in $(seq 100); do
        splitfademographic_analysis/${SAMPLE_NAME}.psmcfa >demographic_analysis/${SAMPLE_NAME}_split.psmcfa
        psmc -b -N25 -t15 -r5 -p "4+25*2+4+6" -odemographic_analysis/bootstrap/${SAMPLE_NAME}_bootstrap${i}.psmcdemographic_analysis/${SAMPLE_NAME}_split.psmcfa
    done
done

# Plot PSMC results
for PSMC_RESULT indemographic_analysis/*.psmc; do
    psmc_plot.pl -u3.75e-8  -g1demographic_analysis/$(basename ${PSMC_RESULT} .psmc) ${PSMC_RESULT}
done

# Step 2: GADMA analysis using dadi engine
python -m gadma --input_vcfsnp_data.vcf --outdirdemographic_analysis/gadma_output --engine dadi --n_clusters 3

# Step 3: Model selection using DIYABC-RF
mkdir -pdemographic_analysis/diyabc_rf
for MODEL in {1..16}; do
    diyabc -pdemographic_analysis/diyabc_rf/model${MODEL} -odemographic_analysis/diyabc_rf/model${MODEL}_results -n 10000 -b 0 -r
done

# Apply Random Forest algorithm to select the best model
Rscript -e "library('DIYABC'); runRFanalysis(path='${OUTPUT_DIR}/diyabc_rf', models=16, output='${OUTPUT_DIR}/best_model.txt')"

# Visualize the optimal model using DemesDraw
python -m demesdraw.utils.plot_model --inputdemographic_analysis/best_model.txt --outputdemographic_analysis/optimal_model_plot.png

# Step 4: Investigate gene flow using Dsuite to compute D-statistics
Dsuite Dtriossnp_data.vcf outgroup_genome.fasta -odemographic_analysis/Dsuite_results
Dsuite Dinvestigatesnp_data.vcf outgroup_genome.fasta -odemographic_analysis/D_f_statistics

# Step 5: Estimate historic gene flow using Migrate-n
migrate-n -idemographic_analysis/migrate_input/migrate_input.nex -odemographic_analysis/migrate_output/migrate_results -b 100000 -m -t 100000

echo "Demographic history and gene flow analysis completed. Results are saved indemographic_analysis."
