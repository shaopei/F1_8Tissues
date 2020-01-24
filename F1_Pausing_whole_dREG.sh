### identify trasncript annotation that have a dREG sites near the 5' 100bp regions
studyBed=dREG
ln -s /workdir/sc2457/F1_Tissues/dREG/Browser/ .
ln -s /workdir/sc2457/F1_Tissues/SingleBaseRunOn/map2ref_1bpbed .

# Use dREG sites with  mat reads >=5 AND pat reads >=5 (not strand specific (NS))
for Head in HT KD SK
do
#  bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}' > ${intermediate_file}
  bedtools coverage   -a <(bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_${studyBed}_5mat5pat_NS_uniq.bed &
done
wait

# label the dREG sites as in dREG100
# output the whole dREG region with plus and minus strand (2 lines)
for Head in HT KD SK
do
  python Fake_Find_span_between_max_read_spots_reportWholedREG.py ${Head}_${studyBed}_5mat5pat_NS_uniq.bed ${Head}_${studyBed}_5mat5pat_NS_uniq_labeled.bed & 
done
wait

# Use dREG sites with  mat reads >=5 AND pat reads >=5 (strand specific)
# one dREG can be plus only, minus only, or both plus and minus with  mat reads >=5 AND pat reads >=5
for Head in HT KD SK
do
#  bedtools coverage -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}' > ${intermediate_file}
  bedtools coverage -s -a <(bedtools coverage -a ${Head}_${studyBed}_5mat5pat_NS_uniq_labeled.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5, $6}'| sort-bed - |uniq > ${Head}_${studyBed}_5mat5pat_uniq.bed &
done
wait


unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered
# generate a smaller SNP file for IGV, SNPs within 100bp of the dREG sites
intersectBed -sorted -u -b <(zcat Browser/*.dREG.peak.score.bed.gz| awk '{OFS="\t"}{print $1, $2-100, $3+100}'|sort-bed -) -a <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }' |sort-bed -) | gzip > SNP_in_dREG.bed.gz &


# identify the abundance of PolII at each position
# strand specific
wait
for Head in HT KD SK
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp.bed &
done
done
wait 

# use python script to generaye input for KS test
for Head in HT KD SK
do
  for allele in mat pat
  do
python Generate_vector_input_for_KStest_NodupPlusMinus.py ${Head}_${allele}_temp.bed ${Head}_${studyBed}_5mat5pat_uniq_${allele}.perBase.bed &
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
  wc -l ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed
done

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
  cat ${Tissue}_${studyBed}_5mat5pat_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($8 <= 0.1){print $0, $1c$2d$3 }'  > ${Tissue}_${studyBed}_5mat5pat_uniq_pValue_fdr0.1.bed
done




