
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
# identify trasncript annotation that have a dREG sites near the 5' end
# keep first 100bp regions
cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN{OFS="\t"}  ($6=="-") {print $1, $3-100, $3, $4, $5, $6, $0}; 
($6=="+") {print $1, $2, $2+100, $4, $5, $6, $0}' > gencode.vM20.annotation_transcript_100bp.bed

for Head in BN HT  SK  SP  LV  GI  ST #KD
do 
intersectBed -a gencode.vM20.annotation_transcript_100bp.bed -b <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz) > ${Head}_gencode.vM20.annotation_transcript.bed
#|awk 'BEGIN{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14}' 
done