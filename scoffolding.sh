fq1=$hicdir/hic_1.fastq.gz
fq2=$hicdir/hic_2.fastq.gz

#Build Index
bwa index  contig.asm.fasta
samtools faidx contig.asm.fasta

## aln
bwa aln -t 6 contig.asm.fasta $fq1 > sample_R1.sai
bwa aln -t 6 contig.asm.fasta $fq2 > sample_R2.sai

bwa sampe contig.asm.fasta sample_R1.sai sample_R2.sai $fq1 $fq2 > sample.bwa_aln.sam

## Filtering 
ALLHiC/scripts/PreprocessSAMs.pl   sample.bwa_aln.sam contig.asm.fasta HINDIII
ALLHiC/scripts/filterBAM_forHiC.pl  sample.bwa_aln.REduced.paired_only.bam sample.clean.sam
samtools view -bt  contig.asm.fasta.fai  sample.clean.sam >  sample.clean.bam

ALLHiC_partition -b sample.clean.bam -r contig.asm.fasta -e HindIII -k 21

## Extract CLM file and counts of restriction sites  #(HindIII: AAGCTT; MboI: GATC)
allhic extract sample.clean.bam contig.asm.fasta --RE AAGCTT

## sort
for K in `seq  1 $k`;
do
    allhic optimize sample.clean.counts_AAGCTT.${k}g${K}.txt sample.clean.clm
done

## final results
ALLHiC_build contig.asm.fasta


#Build Index
samtools faidx  groups.asm.fasta

cut -f 1-2 groups.asm.fasta.fai |grep "sample.clean.counts" > groups.chrn.list

#plot
ALLHiC_plot -b sample.clean.bam -a groups.agp -l groups.chrn.list -m 500k -o groups.hic.heatmap
