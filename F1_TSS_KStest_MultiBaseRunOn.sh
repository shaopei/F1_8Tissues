cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSS_KStest_MultiBaseRunOn

# TSS identify see https://github.com/shaopei/F1_8Tissues/blob/master/F1_TSN_identifyTSS_MultiBaseRunOn.sh
# TSS are strand-specific
ln -s ../identifyTSS_MultiBaseRunOn/*_allReads_TSS.bed .
studyBed=allReads_TSS

# Use TSS sites with  mat reads >=5 AND pat reads >=5 (strand specific) 
# remove TSS in chrX
ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -s -a <(bedtools coverage -sorted -a <(grep -v chrX ${Head}_${studyBed}.bed) -b <(zcat map2ref_1bpbed_map5/${Head}_*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz | sort-bed - --max-mem 10G ) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') \
   -b <(zcat map2ref_1bpbed_map5/${Head}_*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz | sort-bed - --max-mem 10G ) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5, $6}'\
   | sort-bed - --max-mem 10G |uniq > ${Head}_${studyBed}_5mat5pat_uniq.bed &
done
wait

[sc2457@cbsudanko TSS_KStest_MultiBaseRunOn]$ wc -l *bed
   78652 BN_allReads_TSS.bed
   10991 BN_allReads_TSS_5mat5pat_uniq.bed
   41125 GI_allReads_TSS.bed
    7180 GI_allReads_TSS_5mat5pat_uniq.bed
   42521 HT_allReads_TSS.bed
    6747 HT_allReads_TSS_5mat5pat_uniq.bed
   39004 KD_allReads_TSS.bed
    6662 KD_allReads_TSS_5mat5pat_uniq.bed
   91606 LV_allReads_TSS.bed
   17595 LV_allReads_TSS_5mat5pat_uniq.bed
   32363 SK_allReads_TSS.bed
    4568 SK_allReads_TSS_5mat5pat_uniq.bed
   53686 SP_allReads_TSS.bed
    9280 SP_allReads_TSS_5mat5pat_uniq.bed
   34172 ST_allReads_TSS.bed
    5223 ST_allReads_TSS_5mat5pat_uniq.bed
  481375 total



# identify the abundance of PolII at each position
# strand specific
# use all reads (not just allelic reads)
wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz |sort-bed - --max-mem 10G  ) > ${Head}_allReads_temp0.bed &
done

# identify the abundance of PolII at each position
# strand specific
# use ONLY allelic reads
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
  bedtools coverage -sorted -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed_map5/${Head}_*_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed - --max-mem 10G) > ${Head}_${allele}_temp0.bed &
done
done
wait 

# paste the all reads  and alleleic reads count into a table
# only keep base with at least b all reads ($8 >= b)
b=2
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
    paste ${Head}_allReads_temp0.bed ${Head}_${allele}_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"} ($2==$10 && $3==$11 && $8+0 >= b ) {print $0}' |cut -f 9- > ${Head}_${allele}_temp.bed &
  done
done

wait
ln -s ../Generate_vector_input_for_KStest_NodupPlusMinus.py .
ln -s ../KStest_flexible_length.R .
# use python script to generaye input for KS test
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
python Generate_vector_input_for_KStest_NodupPlusMinus.py ${Head}_${allele}_temp.bed ${Head}_${studyBed}_5mat5pat_uniq_${allele}.perBase.bed &
done
done
wait

# get p-value for KS test in R
for Head in BN HT  SK  SP  KD  LV  GI  ST 
do
R --vanilla --slave --args $(pwd) ${Head} ${studyBed}_5mat5pat_uniq < KStest_flexible_length.R &
# output ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed bed6, col7 p.value, col8 fdr
done

# exmaine TSS in IGV
f=0.1
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= f){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_fdr${f}.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= f && $6 =="+"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_fdr${f}_plus.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= f && $6 =="-"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_fdr${f}_minus.bed
done

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


###

### seperate TSS with some kind of allelic difference into two group: 1. change driven by 1 single base 2. Change driven by more than 1 base
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/maxTSN_TSS_TID_combine_analysis_MultiBaseRunOn
# combine asTSS and TSS(as or not as) with as.maxTSN
#${Head}_${studyBed}_5mat5pat_uniq.bed is ${Head}_allReads_TSS_5mat5pat_uniq.bed

for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  # identify TSS_with_as.maxTSN -a TSS -b as.maxTSN
  bedtools intersect -u -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b ${Head}_${studyBed}_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_IGV.bed > ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSN.bed
  # combine TSS_with_as.maxTSN  and asTSS
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSN.bed ${Head}_${studyBed}_5mat5pat_uniq_fdr0.1.bed| awk '{OFS="\t"} {split($4,a,":"); print $1, $2, $3, a[1], $5, $6}' |sort-bed --max-mem 10G - |uniq > ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1.bed
  wc -l ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSN.bed ${Head}_${studyBed}_5mat5pat_uniq_fdr0.1.bed ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1.bed
  # 623 BN_allReads_TSS_5mat5pat_uniq_WithAsMaxTSN.bed
  # 1330 BN_allReads_TSS_5mat5pat_uniq_fdr0.1.bed
  # 1684 BN_allReads_TSS_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1.bed
done

## mask the TSN with biggest difference between the two alleles
# identify the abundance of PolII at each position
# strand specific
# use all reads (not just allelic reads)
ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -d -s -a ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz |sort-bed - --max-mem 10G  ) > ${Head}_allReads_temp0.bed  &
done

# identify the abundance of PolII at each position
# strand specific
# use ONLY allelic reads
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
  bedtools coverage -sorted -d -s -a ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1.bed -b <(zcat map2ref_1bpbed_map5/${Head}_*_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed - --max-mem 10G) > ${Head}_${allele}_temp0.bed &
done
done
wait 

# paste the all reads  and alleleic reads count into a table
# only keep base with at least b all reads ($8 >= b)
b=2
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
    paste ${Head}_allReads_temp0.bed ${Head}_${allele}_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"} ($2==$10 && $3==$11 && $8+0 >= b ) {print $0}' |cut -f 9- > ${Head}_${allele}_temp.bed &
  done
done
wait

# mask the TSN with biggest difference between the two alleles
# for TSS with only 1bp length, it will be removed after the masking
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
paste ${Head}_mat_temp.bed ${Head}_pat_temp.bed  | awk -v b=$b 'BEGIN{OFS="\t"} ($2==$10 && $3==$11 && $7 == $15 ) {print $0}'> ${Head}_mat_pat_temp.bed 
python MaskTSNwithMaxAllelicDifference.py ${Head}_mat_pat_temp.bed ${Head}_mat_pat_temp_masked.bed 
cat ${Head}_mat_pat_temp_masked.bed | cut -f 1-8 > ${Head}_mat_temp_masked.bed
cat ${Head}_mat_pat_temp_masked.bed | cut -f 9-16 > ${Head}_pat_temp_masked.bed
done

wait
ln -s ../Generate_vector_input_for_KStest_NodupPlusMinus.py .
ln -s ../KStest_flexible_length.R .
# use python script to generaye input for KS test
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
python Generate_vector_input_for_KStest_NodupPlusMinus.py ${Head}_${allele}_temp_masked.bed ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_${allele}.perBase.bed &
done
done
wait

# get p-value for KS test in R
# the input bed region that do not have at least 5 pat AND 5 mat reads will be removed in KStest_flexible_length.R
for Head in BN HT  SK  SP  KD  LV  GI  ST 
do
cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_mat.perBase.bed | cut -f 1-6 > ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked.bed 
R --vanilla --slave --args $(pwd) ${Head} ${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked < KStest_flexible_length.R &
# output ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed bed6, col7 p.value, col8 fdr
done

# exmaine TSS in IGV
pValue=0.05
f=0.1
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= 0.1){print $0 }' | wc -l
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v p=$pValue 'BEGIN{OFS="\t"; c=":"; d="-"} ($7+0 <= 0.05){print $0 }' | wc -l
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 >= 0.9){print $0 }' | wc -l
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v p=$pValue 'BEGIN{OFS="\t"; c=":"; d="-"} ($7+0 >0.5){print $0 }' | wc -l

  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} {print $0, $1c$2d$3 }' >  ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue_check.bed 
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} {print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - >  ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue_IGV.bed 
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($6 =="+"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue_IGV_plus.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($6 =="-"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue_IGV_minus.bed
done

# make a table to comapre masked and unmasked p-value and fdr
ln -s ../TSS_KStest_MultiBaseRunOn/*_${studyBed}_5mat5pat_uniq_pValue.bed .
# -a masked -b umasked
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
bedtools intersect -wb -s -a ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_masked_pValue.bed \
-b ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"} {print $9,$10,$11,$12,$13,$14,$7,$8,$15,$16}' \
> ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue.bed 
# output is bed6, col7 maksed pavlue,col 8 masked fdr,col 9 unmaksed pvalue, col 10 unmaked fdr

# as.TSS driven by muliple base (p-value <= cut off) cut off determine by the p-value of the unmaksed KS test (fdr<=0.1)
f=0.1
cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"} ($10+0 <= f){print $0}' > ${Head}.temp2
P_cutoff=`R --vanilla --slave --args ${Head}.temp2 9 < getColMax.R`
cat ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue.bed | awk -v p=$P_cutoff 'BEGIN{OFS="\t"} ($7+0 <= p){print $0}' > ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValueLessThanCutOff${P_cutoff}.bed 


### examine SNPs distribution in asTSS: driven by single base VS driven by multiple base
#bed6, col7 maksed pavlue,col 8 masked fdr,col 9 unmaksed pvalue, col 10 unmaked fdr 
#${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue.bed 

# label the location of maxTSS
 ln -s ../identifyTSS_maxTSNs_MultiBaseRunOn/*_allReads_TSS_maxTSNs.bed .
 bedtools intersect -wo -s -a ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue.bed -b ${Head}_allReads_TSS_maxTSNs.bed > ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue_maxTSNs.bed 
# distribuiton of SNPs in TSS as comparation
 bedtools intersect -wo -s -a <(cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed |awk '{OFS="\t"} {print $1, $2,$3,$4, $5,$6, "NA", "NA", $7,$8 }' ) -b ${Head}_allReads_TSS_maxTSNs.bed > ${Head}_${studyBed}_5mat5pat_uniq_pValue_maxTSNs.bed 
Rscript getSNPsAbundance.R













#### remove below later
### use ${Head}_${studyBed}_5mat5pat_uniq.bed as input, not ${Head}_${studyBed}_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1.bed 
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSS_KStest_MultiBaseRunOn
# mask the TSN with biggest difference between the two alleles
ln -s ../MaskTSNwithMaxAllelicDifference.py .
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
paste ${Head}_mat_temp.bed ${Head}_pat_temp.bed  | awk -v b=$b 'BEGIN{OFS="\t"} ($2==$10 && $3==$11 && $7 == $15 ) {print $0}'> ${Head}_mat_pat_temp.bed 
python MaskTSNwithMaxAllelicDifference.py ${Head}_mat_pat_temp.bed ${Head}_mat_pat_temp_masked.bed 
cat ${Head}_mat_pat_temp_masked.bed | cut -f 1-8 > ${Head}_mat_temp_masked.bed
cat ${Head}_mat_pat_temp_masked.bed | cut -f 9-16 > ${Head}_pat_temp_masked.bed
done

wait
ln -s ../Generate_vector_input_for_KStest_NodupPlusMinus.py .
ln -s ../KStest_flexible_length.R .
# use python script to generaye input for KS test
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
python Generate_vector_input_for_KStest_NodupPlusMinus.py ${Head}_${allele}_temp_masked.bed ${Head}_${studyBed}_5mat5pat_uniq_masked_${allele}.perBase.bed &
done
done
wait

# get p-value for KS test in R
for Head in BN HT  SK  SP  KD  LV  GI  ST 
do
  ln -s  ${Head}_${studyBed}_5mat5pat_uniq.bed ${Head}_${studyBed}_5mat5pat_uniq_masked.bed
R --vanilla --slave --args $(pwd) ${Head} ${studyBed}_5mat5pat_uniq_masked < KStest_flexible_length.R &
# output ${Head}_${studyBed}_5mat5pat_uniq_masked_pValue.bed bed6, col7 p.value, col8 fdr
done

# exmaine TSS in IGV
f=0.1
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_masked_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= f){print $0 }' | wc -l
done
