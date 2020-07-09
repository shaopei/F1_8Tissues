
# merge mat, pat, iden to all reads bigwig track
export mouse_genome=/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz
export mouse_chinfo=/local/storage/data/mm10/mm10.chromInfo

for Head in  BN HT  SK  SP  KD  LV  GI  ST
do
# Convert to bedGraph ... 
j=${Head}_map2ref_1bpbed_map5
# merge mat, pat, iden to all reads bigwig track
bedtools genomecov -bg -i <(zcat ${Head}_*.map2ref.map5.1bp.sorted.bed.gz |LC_COLLATE=C sort -k 1,1) -g ${mouse_chinfo} -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_plus.bedGraph &
bedtools genomecov -bg -i <(zcat ${Head}_*.map2ref.map5.1bp.sorted.bed.gz |LC_COLLATE=C sort -k 1,1) -g ${mouse_chinfo} -strand - |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_minus.bedGraph &
done
wait
# Then to bigWig
for Head in  BN HT  SK  SP  KD  LV  GI  ST
do
  j=${Head}_map2ref_1bpbed_map5
  cat $j\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $j\_minus.inv.bedGraph 
  bedGraphToBigWig $j\_minus.inv.bedGraph ${mouse_chinfo} $j\_minus.bw &
  bedGraphToBigWig $j\_plus.bedGraph ${mouse_chinfo} $j\_plus.bw &
done


# make bigwig tracks from mat, pat, iden seperatley 
for f in *bowtie.gz_AMBremoved_sorted_*.map2ref.map5.1bp.sorted.bed.gz
do 
j=`echo $f |rev |cut -d . -f 4-|rev`
echo $j
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_plus.bedGraph &
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand - |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_minus.bedGraph &
done
wait
for f in *bowtie.gz_AMBremoved_sorted_*.map2ref.map5.1bp.sorted.bed.gz
do 
j=`echo $f |rev |cut -d . -f 4-|rev`
  cat $j\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $j\_minus.inv.bedGraph 
  bedGraphToBigWig $j\_minus.inv.bedGraph ${mouse_chinfo} $j\_minus.bw &
  bedGraphToBigWig $j\_plus.bedGraph ${mouse_chinfo} $j\_plus.bw &
done

for Head in  BN HT  SK  SP  KD  LV  GI  ST
do
  for strand in plus minus
  do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_map2ref_1bpbed_map5_B6_${strand}.bw ${Head}_*B6_all_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp_${strand}.bw  &"
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_map2ref_1bpbed_map5_CAST_${strand}.bw ${Head}_*B6_all_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp_${strand}.bw & "
done
done

