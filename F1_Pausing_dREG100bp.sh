### identify trasncript annotation that have a dREG sites near the 5' 100bp regions
studyBed=dREG

# Use dREG sites with  mat reads >=5 AND pat reads >=5 (not strand specific)
for Head in HT KD SK
do
#  bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}' > ${intermediate_file}
  bedtools coverage -a <(bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_${studyBed}_5mat5pat_uniq.bed &
done
wait

# identify the position with max read count on plus strand and minus strand
# find the mid points between the points with max read count
# generate bed files with plus = mid points + 100, minus = mid points -100 bp 
for Head in HT KD SK
do
  bedtools coverage -d -s -a <(cat ${Head}_${studyBed}_5mat5pat_uniq.bed |awk 'BEGIN{OFS="\t"} {print $0}{print $1, $2, $3, $4, ".", "-"}' ) -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.*at.bowtie.gz_AMBremoved_sorted_*.map2ref.1bp.sorted.bed.gz ) > ${Head}_temp.bed &
done
wait
for Head in HT KD SK
do
  python Find_span_between_max_read_spots.py ${Head}_temp.bed ${Head}_${studyBed}_100bp.bed & 
done
wait

# keep the regions with at least 5 mat and 5 pat reads (strand-specific)
for Head in HT KD SK
do
  bedtools coverage -s -a <(bedtools coverage -s -a ${Head}_${studyBed}_100bp.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_${studyBed}_100bp_5mat5pat_uniq.bed &
done
wait



# identify the abundance of PolII at each position within the bed file generate above
# strand specific
for Head in HT KD SK
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_${studyBed}_100bp_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp.bed &
done
done
wait

# use python script to generaye input for KS test
for Head in HT KD SK
do
  for allele in mat pat
  do
python Generate_vector_input_for_KStest_NodupPlusMinus.py ${Head}_${allele}_temp.bed ${Head}_${studyBed}_100bp_5mat5pat_uniq_${allele}.perBase.bed &
done
done
wait


# get p-value for KS test in R
for Tissue in KD SK HT 
do
R --vanilla --slave --args $(pwd) ${Tissue} ${studyBed}_100bp_5mat5pat_uniq < KStest_flexible_length.R &
done

# HERE!!!
for Head in HT KD SK
do
  echo ${Head}_all.dREG.peak.score.bed.gz
zcat Browser/${Head}_all.dREG.peak.score.bed.gz |wc -l
done
wait
for Tissue in HT KD SK
do
  cat ${Tissue}_${studyBed}_100bp_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($7 <= 0.05){print $0, $1c$2d$3 }' |wc -l
done
for Tissue in HT KD SK
do
  cat ${Tissue}_${studyBed}_100bp_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0, $1c$2d$3 }' |wc -l
done

for Tissue in HT KD SK
do
  echo ${Tissue}_${studyBed}_100bp_5mat5pat_uniq_pValue.bed 
  cat ${Tissue}_${studyBed}_100bp_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0, $1c$2d$3 }' > ${Tissue}_${studyBed}_100bp_5mat5pat_uniq_pValue_fdr0.1.bed
done




