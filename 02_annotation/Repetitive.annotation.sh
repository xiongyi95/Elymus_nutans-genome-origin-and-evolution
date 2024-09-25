#!/bin/bash

# Define paths and parameters
GENOME="path/to/genome.fasta"  # Path to the genome assembly file
REPBASE_LIB="path/to/repbase.lib"  # Path to Repbase library file for RepeatMasker
OUTPUT_DIR="repeats_annotation"  # Output directory for the results
LTR_DB="Enutans_LTR_DB.fasta"  # Output file for the custom LTR-RT database

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Evidence-based annotation using RepeatMasker with the Repbase database
RepeatMasker -pa 8 -lib ${REPBASE_LIB} -dir ${OUTPUT_DIR}/repeatmasker_output ${GENOME}

# Step 2: De novo LTR-RT database construction
# Step 2.1: Run LTRharvest to identify LTR elements
gt suffixerator -db ${GENOME} -indexname ${OUTPUT_DIR}/genome_index -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ${OUTPUT_DIR}/genome_index -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -similar 85 -vic 10 -seed 20 \
    -seqids yes -out ${OUTPUT_DIR}/ltrharvest_output.gff3 -gff3 yes

# Step 2.2: Run LTR_Finder to find LTR elements
LTR_Finder ${GENOME} -s 85 -L 7000 -M 4 -f 10 -g ${OUTPUT_DIR}/ltrfinder_output.gff3

# Step 2.3: Consolidate results using LTR_retriever
LTR_retriever -genome ${GENOME} -inharvest ${OUTPUT_DIR}/ltrharvest_output.gff3 -infinder ${OUTPUT_DIR}/ltrfinder_output.gff3 -threads 8 -o ${OUTPUT_DIR}/ltr_retriever_output

# Step 3: Combine LTRharvest, LTR_Finder, and RepeatMasker results to create a non-redundant LTR-RT database
cat ${OUTPUT_DIR}/ltr_retriever_output/*.fasta > ${OUTPUT_DIR}/${LTR_DB}

# Step 4: Run RepeatMasker with the custom LTR-RT database
RepeatMasker -pa 8 -lib ${OUTPUT_DIR}/${LTR_DB} -dir ${OUTPUT_DIR}/ltr_repeatmasker_output ${GENOME}

# Step 5: Combine evidence-based and de novo predictions into the final repetitive element annotation
cat ${OUTPUT_DIR}/repeatmasker_output/*.out ${OUTPUT_DIR}/ltr_repeatmasker_output/*.out > ${OUTPUT_DIR}/final_repeats_annotation.out

# Clean up intermediate files (optional)
# rm -r ${OUTPUT_DIR}/repeatmasker_output
# rm -r ${OUTPUT_DIR}/ltr_repeatmasker_output
# rm -r ${OUTPUT_DIR}/genome_index*

echo "Repetitive elements annotation completed. Final results are in ${OUTPUT_DIR}/final_repeats_annotation.out"
