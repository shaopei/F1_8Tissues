cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TID.dREG_KStest_MultiBaseRunOn

studyBed=dREG
ln -s /workdir/sc2457/F1_Tissues/dREG/Browser/ .
ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .
ln -s ../Fake_Find_span_between_max_read_spots_reportWholedREG.py .
ln -s ../Generate_vector_input_for_KStest_NodupPlusMinus.py .
ln -s ../KStest_flexible_length.R .

# Use dREG sites with  mat reads >=5 AND pat reads >=5 (not strand specific (NS))
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage   -a <(bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"}') -b <(zcat map2ref_1bpbed_map5/${Head}_*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed_map5/${Head}_*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_${studyBed}_5mat5pat_NS_uniq.bed &
done
wait

[sc2457@cbsudanko TID.dREG_KStest_MultiBaseRunOn]$ wc *bed -l
  21420 BN_dREG_5mat5pat_NS_uniq.bed
  16950 GI_dREG_5mat5pat_NS_uniq.bed
  11426 HT_dREG_5mat5pat_NS_uniq.bed
  16600 KD_dREG_5mat5pat_NS_uniq.bed
  26743 LV_dREG_5mat5pat_NS_uniq.bed
  10271 SK_dREG_5mat5pat_NS_uniq.bed
  16807 SP_dREG_5mat5pat_NS_uniq.bed
   9311 ST_dREG_5mat5pat_NS_uniq.bed
 129528 total

# label the dREG sites as in dREG100
# output the whole dREG region with plus and minus strand (2 lines)
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  python Fake_Find_span_between_max_read_spots_reportWholedREG.py ${Head}_${studyBed}_5mat5pat_NS_uniq.bed ${Head}_${studyBed}_5mat5pat_NS_uniq_labeled.bed & 
done
wait

# Use dREG sites with  mat reads >=5 AND pat reads >=5 (strand specific) (above was non-strand specific, here furthur restrict 5 mat reads AND 5 pat reads per strand)
# one dREG can be plus only, minus only, or both plus and minus with  mat reads >=5 AND pat reads >=5
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -s -a <(bedtools coverage -a ${Head}_${studyBed}_5mat5pat_NS_uniq_labeled.bed -b <(zcat map2ref_1bpbed_map5/${Head}_*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed_map5/${Head}_*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5, $6}'| sort-bed - |uniq > ${Head}_${studyBed}_5mat5pat_uniq.bed &
done
wait

[sc2457@cbsudanko TID.dREG_KStest_MultiBaseRunOn]$ wc -l *bed
   21420 BN_dREG_5mat5pat_NS_uniq.bed # keep those with mat reads >=5 AND pat reads >=5 (NON-strand specific)
   42840 BN_dREG_5mat5pat_NS_uniq_labeled.bed # twice(+ and - strand) of _NS_uniq.bed
   25437 BN_dREG_5mat5pat_uniq.bed  # keep those with mat reads >=5 AND pat reads >=5 (strand specific)
   16950 GI_dREG_5mat5pat_NS_uniq.bed
   33900 GI_dREG_5mat5pat_NS_uniq_labeled.bed
   19840 GI_dREG_5mat5pat_uniq.bed
   11426 HT_dREG_5mat5pat_NS_uniq.bed
   22852 HT_dREG_5mat5pat_NS_uniq_labeled.bed
   13359 HT_dREG_5mat5pat_uniq.bed
   16600 KD_dREG_5mat5pat_NS_uniq.bed
   33200 KD_dREG_5mat5pat_NS_uniq_labeled.bed
   19333 KD_dREG_5mat5pat_uniq.bed
   26743 LV_dREG_5mat5pat_NS_uniq.bed
   53486 LV_dREG_5mat5pat_NS_uniq_labeled.bed
   33956 LV_dREG_5mat5pat_uniq.bed
   10271 SK_dREG_5mat5pat_NS_uniq.bed
   20542 SK_dREG_5mat5pat_NS_uniq_labeled.bed
   12160 SK_dREG_5mat5pat_uniq.bed
   16807 SP_dREG_5mat5pat_NS_uniq.bed
   33614 SP_dREG_5mat5pat_NS_uniq_labeled.bed
   20516 SP_dREG_5mat5pat_uniq.bed
    9311 ST_dREG_5mat5pat_NS_uniq.bed
   18622 ST_dREG_5mat5pat_NS_uniq_labeled.bed
   10731 ST_dREG_5mat5pat_uniq.bed
  543916 total


unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered
# generate a smaller SNP file for IGV, SNPs within 100bp of the dREG sites
intersectBed -sorted -u -b <(zcat Browser/*.dREG.peak.score.bed.gz| awk '{OFS="\t"}{print $1, $2-100, $3+100}'|sort-bed -) -a <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }' |sort-bed -) | gzip > SNP_in_dREG.bed.gz &

# identify the abundance of PolII at each position
# strand specific
# use all reads (not just allelic reads)
wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz ) > ${Head}_allReads_temp0.bed &
done

# identify the abundance of PolII at each position
# strand specific
# use ONLY allelic reads

for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed_map5/${Head}_*_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp0.bed &
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
    paste ${Head}_allReads_temp0.bed ${Head}_${allele}_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"} ($4==$12 && $8 >= b ) {print $0}' |cut -f 9- > ${Head}_${allele}_temp.bed &
  done
done

wait
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
done

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



for Tissue in HT KD SK
do
  echo ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed 
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0, $1c$2d$3 }'  > ${Tissue}_${studyBed}_5mat5pat_uniq_pValue_fdr0.1.txt
done

for Tissue in HT KD SK
do
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Tissue}_${studyBed}_5mat5pat_uniq_pValue_fdr0.1_TSS.bed
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($6 =="+"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Tissue}_${studyBed}_5mat5pat_uniq_pValue_TSS_plus.bed
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($6 =="-"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Tissue}_${studyBed}_5mat5pat_uniq_pValue_TSS_minus.bed
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0}' |sort-bed - > ${Tissue}_${studyBed}_5mat5pat_uniq_pValue_fdr0.1.bed
done


for Tissue in HT KD SK
do
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 > 0.9){print $0}' > ${Tissue}_${studyBed}_5mat5pat_uniq_pValue_fdr0.9.bed &
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


