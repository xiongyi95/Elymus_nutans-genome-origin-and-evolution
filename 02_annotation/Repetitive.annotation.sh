#!/bin/bash

# Step 1: Evidence-based annotation using RepeatMasker with the Repbase database
RepeatMasker -pa 8 -lib repbase.lib -dir ./repeatmasker_output genome.fasta

# Step 2: De novo LTR-RT database construction
# Step 2.1: Run LTRharvest to identify LTR elements
gt suffixerator -db genome.fasta -indexname ./genome_index -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ./genome_index -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -similar 85 -vic 10 -seed 20 \
    -seqids yes -out ./ltrharvest_output.gff3 -gff3 yes

# Step 2.2: Run LTR_Finder to find LTR elements
LTR_Finder genome.fasta -s 85 -L 7000 -M 4 -f 10 -g ./ltrfinder_output.gff3

# Step 2.3: Consolidate results using LTR_retriever
LTR_retriever -genome genome.fasta -inharvest ./ltrharvest_output.gff3 -infinder ./ltrfinder_output.gff3 -threads 8 -o ./ltr_retriever_output

# Step 3: Combine LTRharvest, LTR_Finder, and RepeatMasker results to create a non-redundant LTR-RT database
cat ./ltr_retriever_output/*.fasta > ./Enutans_LTR_DB.fasta

# Step 4: Run RepeatMasker with the custom LTR-RT database
RepeatMasker -pa 8 -lib ./Enutans_LTR_DB.fasta -dir ./ltr_repeatmasker_output genome.fasta

# Step 5: Combine evidence-based and de novo predictions into the final repetitive element annotation
cat ./repeatmasker_output/*.out ./ltr_repeatmasker_output/*.out > ./final_repeats_annotation.out

# Clean up intermediate files (optional)
# rm -r ./repeatmasker_output
# rm -r ./ltr_repeatmasker_output
# rm -r ./genome_index*

echo "Repetitive elements annotation completed. Final results are in ./final_repeats_annotation.out"
