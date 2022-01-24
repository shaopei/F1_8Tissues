
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


#make mat, pat bigwig for each sample for heatmap
# BN
cd /workdir/sc2457/F1_Tissues/second_batch/map2ref_1bpbed_map5
ln -s ../allelicbias-PersonalGenome_P.CAST_M.B6-BN_*R1/BN_*_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz .
ln -s ../allelicbias-PersonalGenome_P.CAST_M.B6-BN_*R1/BN_*_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz .

export mouse_genome=/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz
cat /local/storage/data/mm10/mm10.chromInfo | awk'{OFS="\t"} {print $1, $2}' > mouse_chinfo.txt

#for Head in  BN HT  SK  SP  KD  LV  GI  ST
#do
# Convert to bedGraph ... 
Head=BN

for s in  MB6_F MB6_G PB6_B PB6_C PB6_D PB6_E
do 
  for p in mat pat 
  do 
    bedtools genomecov -bg -i <(zcat ${Head}_${s}_R1.${p}.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz | awk '{OFS="\t"}($6=="+") {print $1, $2, $2+1, "o", $5, $6}; ($6=="-") {print $1, $3-1, $3, "o", $5, $6} '|LC_COLLATE=C sort -k 1,1) \
    -g mouse_chinfo.txt -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${Head}_${s}_R1.${p}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted_plus.bedGraph &
    bedtools genomecov -bg -i <(zcat ${Head}_${s}_R1.${p}.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz | awk '{OFS="\t"}($6=="+") {print $1, $2, $2+1, "o", $5, $6}; ($6=="-") {print $1, $3-1, $3, "o", $5, $6} '|LC_COLLATE=C sort -k 1,1) \
    -g mouse_chinfo.txt -strand - |LC_COLLATE=C sort -k1,1 -k2,2n > ${Head}_${s}_R1.${p}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted_minus.bedGraph &
done
wait

for f in *bowtie.gz_AMBremoved_sorted_*.map2ref.sorted.bed.gz 
do 
h=`echo $f |rev |cut -d . -f 4-|rev`
j=${h}.map5.1bp.sorted
  cat $j\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $j\_minus.inv.bedGraph 
  bedGraphToBigWig $j\_minus.inv.bedGraph ${mouse_chinfo} $j\_minus.bw &
  bedGraphToBigWig $j\_plus.bedGraph ${mouse_chinfo} $j\_plus.bw &
done


###

# single base runon
# make bigwig tracks from mat, pat, iden seperatley 
for f in *bowtie.gz_AMBremoved_sorted_*.map2ref.map5.1bp.sorted.bed.gz
do 
j=`echo $f |rev |cut -d . -f 4-|rev`
echo $j
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_plus.bedGraph &
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand - |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_minus.inv.bedGraph &
done
wait
for f in *bowtie.gz_AMBremoved_sorted_*.map2ref.map5.1bp.sorted.bed.gz
do 
j=`echo $f |rev |cut -d . -f 4-|rev`
  cat $j\_minus.inv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $j\_minus.bedGraph 
  bedGraphToBigWig $j\_minus.bedGraph ${mouse_chinfo} $j\_minus.bw &
  bedGraphToBigWig $j\_plus.bedGraph ${mouse_chinfo} $j\_plus.bw &
done

for Head in  HT  SK  KD
do
  for strand in plus minus
  do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_PB6_F5N6_map2ref_1bpbed_map5_B6_${strand}.bw map2ref_1bpbed_map5/${Head}_PB6_*_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp_${strand}.bw  &"
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_PB6_F5N6_map2ref_1bpbed_map5_CAST_${strand}.bw map2ref_1bpbed_map5/${Head}_PB6_*_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp_${strand}.bw & "
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_PB6_F5N6_map2ref_1bpbed_map5_ALL_${strand}.bw map2ref_1bpbed_map5/${Head}_PB6_*_R1.*at.bowtie.gz_AMBremoved_sorted_*.map2ref.map5.1bp_${strand}.bw & "
done                                                                                                  
done




# map3
for f in *bowtie.gz_AMBremoved_sorted_*.map2ref.1bp.sorted.bed.gz
do 
j=`echo $f |rev |cut -d . -f 4-|rev`
echo $j
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_plus.bedGraph &
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand - |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_minus.inv.bedGraph &
done
wait
for f in *bowtie.gz_AMBremoved_sorted_*.map2ref.1bp.sorted.bed.gz
do 
j=`echo $f |rev |cut -d . -f 4-|rev`
  echo $j
  cat $j\_minus.inv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $j\_minus.bedGraph 
  bedGraphToBigWig $j\_minus.bedGraph ${mouse_chinfo} $j\_minus.bw &
  bedGraphToBigWig $j\_plus.bedGraph ${mouse_chinfo} $j\_plus.bw &
done

for Head in  HT  SK  KD
do
  for strand in plus minus
  do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_PB6_F5N6_map2ref_1bpbed_map3_B6_${strand}.bw map2ref_1bpbed/${Head}_PB6_*_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp_${strand}.bw  &"
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_PB6_F5N6_map2ref_1bpbed_map3_CAST_${strand}.bw map2ref_1bpbed/${Head}_PB6_*_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp_${strand}.bw & "
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Head}_PB6_F5N6_map2ref_1bpbed_map3_ALL_${strand}.bw map2ref_1bpbed/${Head}_PB6_*_R1.*at.bowtie.gz_AMBremoved_sorted_*.map2ref.1bp_${strand}.bw & "
done                                                                                                  
done


# map3 from shared allelic maxTSNs
for Head in SK KD #HT
do
for allele in mat pat
  do 
f=${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads.bed
j=${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_plus.bedGraph &
bedtools genomecov -bg -i $f -g ${mouse_chinfo} -strand - |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_minus.inv.bedGraph &
done
done

wait
for Head in SK KD #HT
do
for allele in mat pat
  do 
f=${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads.bed
j=${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads
  cat $j\_minus.inv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $j\_minus.bedGraph 
  bedGraphToBigWig $j\_minus.bedGraph ${mouse_chinfo} $j\_minus.bw &
  bedGraphToBigWig $j\_plus.bedGraph ${mouse_chinfo} $j\_plus.bw &
done
done
