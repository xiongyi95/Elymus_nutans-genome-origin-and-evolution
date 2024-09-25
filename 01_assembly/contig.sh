#01 pacbio_hifi下机bam转fasta.gz
bam2fasta hifi.bam -c 9 -o sample -j 20 &

#02 hifisam组装基因组
 
hifiasm -o contig.asm   \
        -t 72  \
        ./sample.fasta
        
awk '/^S/{print ">"$2; print $3}' contig.asm.bp.p_ctg.gfa > contig.asm.fasta
