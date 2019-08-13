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
  for allele in mat pat
  do
  bedtools coverage -d -s -a <(cat ${Head}_${studyBed}_5mat5pat_uniq.bed |awk 'BEGIN{OFS="\t"} {print $0}{print $1, $2, $3, $4, ".", "-"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp.bed
  python Find_span_between_max_read_spots.py ${Head}_${allele}_temp.bed ${Head}_${studyBed}_200bp.bed & #${Head}_${studyBed}_200bp_5mat5pat_uniq.bed &
  done
done
wait

# keep the regions with at least 5 mat and 5 pat reads (strand-specific)
for Head in HT KD SK
do
  bedtools coverage -s -a <(bedtools coverage -s -a ${Head}_${studyBed}_200bp.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_${studyBed}_200bp_5mat5pat_uniq.bed &
done
wait



# identify the abundance of PolII at each position within the bed file generate above
# strand specific
for Head in HT KD SK
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_${studyBed}_200bp_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp.bed
    #|cut -f 8| paste - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > ${Head}_${studyBed}_wSNP_5mat5pat+reads_uniq_${allele}.perBase.bed &
  cat ${Head}_${allele}_temp.bed  |cut -f 8| paste - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > ${Head}_${studyBed}_200bp_5mat5pat_uniq_${allele}.perBase.bed &
done
done
wait

# get p-value for KS test in R
for Tissue in KD SK HT 
do
R --vanilla --slave --args $(pwd) ${Tissue} ${studyBed}_200bp_5mat5pat_uniq < KStest.R &
done
wait
# HERE!!!
for Tissue in HT KD SK
do
  echo ${Tissue}_gencode.vM20.annotation_transcript_100bp_5mat5pat_uniq_pValue.bed 
  cat ${Tissue}_gencode.vM20.annotation_transcript_100bp_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($7 <= 0.15){print $0, $1c$2d$3 }'
done


# how many of them are significant? (p<0.05)
for Tissue in HT KD SK
do
  cat ${Tissue}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue.bed |awk 'BEGIN {OFS="\t"} ($7<=0.05) {print $0}' >  ${Tissue}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue0.05.bed &
done

# see the allele-specific polII distribution at genome browser
CHINFO=/local/storage/data/mm10/mm10.chromInfo
for f in *at.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz
do j=`echo $f|rev | cut -d . -f 3-|rev`
   echo $j
   ## Convert to bedGraph ...
   bedtools genomecov -bg -i ${j}.bed.gz -g ${CHINFO} -strand + | sort-bed - > $j\_plus.bedGraph &
   bedtools genomecov -bg -i ${j}.bed.gz -g ${CHINFO} -strand - | sort-bed - > $j\_minus.noinv.bedGraph &
   wait
  ## Invert minus strand.
   cat $j\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > $j\_minus.bedGraph ## Invert read counts on the minus strand.
   ## Then to bigWig (nomalized and non-nomrmalized ones)
   bedGraphToBigWig $j\_plus.bedGraph ${CHINFO} $j\_plus.bw  &
   bedGraphToBigWig $j\_minus.bedGraph ${CHINFO} $j\_minus.bw  &
done
