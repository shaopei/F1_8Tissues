
PREFIX=map2ref_bed/BN_MB6_all_R1
MAT_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
PAT_READ_BED=${PREFIX}.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz

bedtools intersect -a <(zcat BN_all.dREG.peak.score.bed.gz) -b <(zcat ${MAT_READ_BED})

map2ref_bed/BN_MB6_all_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz

for Head in BN HT  SK  SP  LG  LV  GI  ST
do 
echo $Head
bed_dir=map2ref_bed
  for f in ${bed_dir}/${Head}_MB6_*_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz  #each samples from MB6 of the same tissue
    do PREFIX=`echo $f|cut -d . -f 1`
    echo $PREFIX
    MAT_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
    PAT_READ_BED=${PREFIX}.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
    #IDENTICAL_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz

    P=`echo $PREFIX| cut -d / -f 2`
    ln -s HMM_bed/combined_cross/T8_HMM_plus.bed ${P}_HMM_plus.bed
    ln -s HMM_bed/combined_cross/T8_HMM_minus.bed ${P}_HMM_minus.bed
    BinomialTest ${P}_HMM_plus.bed ${P}_HMM_plus ${MAT_READ_BED} ${PAT_READ_BED} &
    BinomialTest ${P}_HMM_minus.bed ${P}_HMM_minus ${MAT_READ_BED} ${PAT_READ_BED} &
    wait
  done

### Ignore above , start from here!!! ####
# identify trasncript annotation that have a dREG sites near the 5' 100bp regions
cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN{OFS="\t"}  ($6=="-") {print $1, $3-100, $3, $4, $5, $6, $0}; 
($6=="+") {print $1, $2, $2+100, $4, $5, $6, $0}' > gencode.vM20.annotation_transcript_100bp.bed

for Head in BN HT  SK  SP  LV  GI  ST KD
do 
intersectBed -a gencode.vM20.annotation_transcript_100bp.bed -b <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz) > ${Head}_gencode.vM20.annotation_transcript_100bp.bed
#|awk 'BEGIN{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14}' 
done

# keep the first 30bp
for Head in BN HT  SK  SP  LV  GI  ST KD
do 
cat ${Head}_gencode.vM20.annotation_transcript_100bp.bed | awk 'BEGIN{OFS="\t"}  ($6=="-") {print $1, $3-30, $3, $4, $5, $6, $0}; 
($6=="+") {print $1, $2, $2+30, $4, $5, $6, $0}'| LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${Head}_gencode.vM20.annotation_transcript_30bp.bed 
done


# Keep TRX with SNPs in the first 30bp
unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered

for Head in BN HT  SK  SP  LV  GI  ST KD
do 
intersectBed -sorted -u -a ${Head}_gencode.vM20.annotation_transcript_30bp.bed -b <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }') > ${Head}_gencode.vM20.annotation_transcript_30bp_wSNP.bed &
#|awk 'BEGIN{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14}' 
done

# Keep TRX with more than 10 reads (sum mat/pat F5/F6) (strand specific)
for Head in HT KD SK
do
  #rm ${Head}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads.bed
  bedtools coverage -s -a <(cat ${Head}_gencode.vM20.annotation_transcript_30bp_wSNP.bed | cut -f 1-6) -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.*at.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz | awk 'BEGIN{OFS="\t"} ($6=="+"){print $1, $2, $3, $4, $5, "-"} ($6=="-"){print $1, $2, $3, $4, $5, "+"}') | awk 'BEGIN{OFS="\t"} ($7 >=10){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq.bed &
done


# the abundance of PolII at each position within the bed file
for Head in HT KD SK
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz | awk 'BEGIN{OFS="\t"} ($6=="+"){print $1, $2, $3, $4, $5, "-"} ($6=="-"){print $1, $2, $3, $4, $5, "+"}') |cut -f 8| paste - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  > ${Head}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_${allele}.perBase.bed &
done
done

# get p-value for KS test in R
for Tissue in HT KD SK
do
R --vanilla --slave --args $(pwd) ${Tissue} < KStest.R &
done

# how many of them are significant? (p<0.05)
for Tissue in HT KD SK
do
  cat ${Tissue}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue.bed |awk 'BEGIN {OFS="\t"} ($7<=0.05) {print $0}' >  ${Tissue}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue0.05.bed &
done




