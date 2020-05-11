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


###



#HERE
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  #echo ${Head}_all.dREG.peak.score.bed.gz
zcat Browser/${Head}_all.dREG.peak.score.bed.gz |wc -l
done
wait

for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  # only use bed regions with at least 5 VALID mat AND 5 pat VALID reads
  wc -l ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed 
done


for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0, $1c$2d$3 }' |wc -l
done


for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 > 0.9){print $0, $1c$2d$3 }' |wc -l
done



for Tissue in BN HT  SK  SP  KD  LV  GI  ST
do
  echo ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed 
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0, $1c$2d$3 }'  > ${Head}_${studyBed}_5mat5pat_uniq_pValue_fdr0.1.txt
done

for Tissue in HT KD SK
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_pValue_fdr0.1_TSS.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($6 =="+"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_pValue_TSS_plus.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($6 =="-"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_pValue_TSS_minus.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0}' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_pValue_fdr0.1.bed
done


for Tissue in HT KD SK
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 > 0.9){print $0}' > ${Head}_${studyBed}_5mat5pat_uniq_pValue_fdr0.9.bed &
done


# get the fasta from dREG sites
# fidn the maxTSN


#HERE

# merge bigwig files from mouse F5 and F6
#cd map2ref_1bpbed
#KD_PB6_F5_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_plus.bw  KD_PB6_F6_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_plus.bw
#KD_PB6_F5_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_plus.bw  KD_PB6_F6_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_plus.bw


 
export mouse_genome=/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz
export mouse_chinfo=/local/storage/data/mm10/mm10.chromInfo
for Tissue in HT KD SK
do
  for parent in mat pat
  do
    for strand in plus minus
    do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Tissue}.${parent}.map2ref.1bp_${strand}.bw map2ref_1bpbed_map5/${Tissue}_PB6_*_dedup_R1.${parent}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_${strand}.bw"
done
done
done

for Tissue in HT KD SK
do
  for parent in mat pat
  do
    for strand in plus minus
    do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Tissue}.${parent}.map5.map2ref.1bp_${strand}.bw map2ref_1bpbed_map5/${Tissue}_PB6_*_dedup_R1.${parent}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted_1bp_${strand}.bw" 
done
done
done

for Tissue in HT KD SK
do
  for parent in mat
  do
    for strand in plus minus
    do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Tissue}.identical.map5.map2ref.1bp_${strand}.bw map2ref_1bpbed_map5/${Tissue}_PB6_*_dedup_R1.${parent}.bowtie.gz_AMBremoved_sorted_identical.map2ref.map5.1bp.sorted_1bp_${strand}.bw" 
done
done
done


#Myproseq2.0Output/KD_PB6_F5_dedup_QC_end_plus.bw  Myproseq2.0Output/KD_PB6_F6_dedup_QC_end_plus.bw
for Tissue in HT KD SK
do
  for strand in plus minus
  do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Tissue}_PB6_F5N6_dedup_QC_end_${strand}.bw Myproseq2.0Output-map3/${Tissue}_PB6_F5_dedup_QC_end_${strand}.bw  Myproseq2.0Output-map3/${Tissue}_PB6_F6_dedup_QC_end_${strand}.bw"
done
done

for Tissue in HT KD SK
do
  for strand in plus minus
  do
    echo "bash mergeBigWigs.bsh --chrom-info=${mouse_chinfo} ${Tissue}_PB6_F5N6_dedup_QC_end_map5_${strand}.bw Myproseq2.0Output-map5/${Tissue}_PB6_F5_dedup_QC_end.sort_${strand}.bw  Myproseq2.0Output-map5/${Tissue}_PB6_F6_dedup_QC_end.sort_${strand}.bw"
done
done

### furthur analysis in R
# make heatmap of proseq reads in short and long pause using 
# SingleRunOn_Pause_analysis.R

### make heatmap of SNPs locations
# make SNPs location bigwig files
ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered .
# generate a smaller SNP file for IGV, SNPs within 100bp of the dREG sites
cat P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6, ".", "+"}' |sort-bed - | bgzip > P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered.bed.gz &
CHINFO=/local/storage/data/mm10/mm10.chromInfo
j=P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered
# Convert to bedGraph ... 
bedtools genomecov -bg -i ${j}.bed.gz -g ${CHINFO} -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_plus.bedGraph 
# Then to bigWig
bedGraphToBigWig $j\_plus.bedGraph ${CHINFO} $j\_plus.bw &


