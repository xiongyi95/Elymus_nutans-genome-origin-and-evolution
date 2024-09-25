#01 pacbio_hifi下机bam转fasta.gz
bam2fasta m84114_231124_090346_s2.hifi_reads.bc2013.bam -c 9 -o sample.2 -j 20 &

#02 hifisam组装基因组
 
hifiasm -o contig.asm   \
        -t 72  \
        ./merged.fasta
        
awk '/^S/{print ">"$2; print $3}' contig.asm.bp.p_ctg.gfa > contig.asm.fasta



