### identify trasncript annotation that have a dREG sites near the 5' 100bp regions
studyBed=dREG

# Use dREG sites with  mat reads >=5 AND pat reads >=5 (strand specific)
for Head in HT KD SK
do
#  bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}' > ${intermediate_file}
  bedtools coverage -s -a <(bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5}'| sort-bed - |uniq |awk 'BEGIN{OFS="\t"} {print $0, "+"} {print $0, "-"}'> ${Head}_${studyBed}_5mat5pat_uniq.bed &
done
wait

# identify the abundance of PolII at each position
# strand specific
for Head in HT KD SK
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp.bed &
done
done
wait

# here add python
for Head in HT KD SK
do
  for allele in mat pat
  do
python Generate_vector_input_for_KStest.py ${Head}_${allele}_temp.bed ${Head}_${studyBed}_5mat5pat_uniq_${allele}.perBase.bed &
done
done
wait

# get p-value for KS test in R
for Tissue in KD SK HT 
do
R --vanilla --slave --args $(pwd) ${Tissue} ${studyBed}_5mat5pat_uniq < KStest_flexible_length.R &
done

# HERE!!!
for Head in HT KD SK
do
  #echo ${Head}_all.dREG.peak.score.bed.gz
zcat Browser/${Head}_all.dREG.peak.score.bed.gz |wc -l
done
wait
for Tissue in HT KD SK
do
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($7 <= 0.05){print $0, $1c$2d$3 }' |wc -l
done
for Tissue in HT KD SK
do
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0, $1c$2d$3 }' |wc -l
done

for Tissue in HT KD SK
do
  echo ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed 
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($7 <= 0.15){print $0, $1c$2d$3 }' 
done




