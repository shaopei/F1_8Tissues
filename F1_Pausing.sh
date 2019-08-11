### identify trasncript annotation that have a dREG sites near the 5' 100bp regions
studyBed=gencode.vM20.annotation_transcript_100bp

# use first 100bp
cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN{OFS="\t"}  ($6=="-") {print $1, $3-100, $3, $4, $5, $6, $0}; 
($6=="+") {print $1, $2, $2+100, $4, $5, $6, $0}' > ${studyBed}.bed


# that intersect with dREG sites
for Head in BN HT  SK  SP  LV  GI  ST KD
do 
intersectBed -wa -a gencode.vM20.annotation_transcript_100bp.bed -b <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz) | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${Head}_${studyBed}.bed
done



unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered
# remove this limit # Keep TRX with SNPs in the first 30bp 
#for Head in BN HT  SK  SP  LV  GI  ST KD
#do 
#intersectBed -sorted -u -a ${Head}_${studyBed}.bed -b <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }') > ${Head}_${studyBed}_wSNP.bed &
#done

# generate a smaller SNP file for IGV, SNPs within 1000bp of the transcript annotation
intersectBed -sorted -u -b <(cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN{OFS="\t"}  ($1 != "chrM"){print $1, $2-1000, $3+1000}' |sort-bed -) -a <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }' |sort-bed -) | gzip > SNP_in_gencode.vM20.annotation_transcript.bed.gz &


# Instead of Keeping TRX with more than 10 reads (sum mat/pat F5/F6) (strand specific)
# do this: new requiements: have mat reads >=5 AND pat reads >=5
for Head in HT KD SK
do
  #intermediate_file=${Head}_${studyBed}_5mat.bed
  #bedtools coverage -s -a <(cat ${Head}_${studyBed}.bed | cut -f 1-6) -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}' > ${intermediate_file}
  bedtools coverage -s -a <(bedtools coverage -s -a <(cat ${Head}_${studyBed}.bed | cut -f 1-6) -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_${studyBed}_5mat5pat_uniq.bed &
done


# the abundance of PolII at each position within the bed file
# sterand specific
# HERE!!!
for Head in HT KD SK
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp.bed  #|cut -f 8| paste - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > ${Head}_${studyBed}_wSNP_5mat5pat+reads_uniq_${allele}.perBase.bed &
  cat ${Head}_${allele}_temp.bed  |cut -f 8| paste - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > ${Head}_${studyBed}_5mat5pat_uniq_${allele}.perBase.bed &
done
done

# get p-value for KS test in R
for Tissue in KD SK HT 
do
R --vanilla --slave --args $(pwd) ${Tissue} < KStest.R &
done

for Tissue in HT KD SK
do
  echo ${Tissue}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue.bed 
  cat ${Tissue}_gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue.bed | awk 'BEGIN{OFS="\t"; c=":"; d="-"} ($7 <= 0.05){print $0, $1c$2d$3 }' 
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
