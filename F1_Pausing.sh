
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

studyBed=gencode.vM20.annotation_transcript_100bp

for Head in BN HT  SK  SP  LV  GI  ST KD
do 
intersectBed -wa -a gencode.vM20.annotation_transcript_100bp.bed -b <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz) | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${Head}_${studyBed}.bed
#|awk 'BEGIN{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14}' 
done


# keep the first 100bp instead of 30bp
#for Head in BN HT  SK  SP  LV  GI  ST KD
#do 
#cat ${Head}_gencode.vM20.annotation_transcript_100bp.bed | awk 'BEGIN{OFS="\t"}  ($6=="-") {print $1, $3-30, $3, $4, $5, $6, $0}; 
#($6=="+") {print $1, $2, $2+30, $4, $5, $6, $0}'| LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${Head}_gencode.vM20.annotation_transcript_30bp.bed 
#done


# remove this limit # Keep TRX with SNPs in the first 30bp 
unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered

#for Head in BN HT  SK  SP  LV  GI  ST KD
#do 
#intersectBed -sorted -u -a ${Head}_${studyBed}.bed -b <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }') > ${Head}_${studyBed}_wSNP.bed &
#done

# generate a smaller SNP file for IGV
intersectBed -sorted -u -b <(cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN{OFS="\t"}  ($1 != "chrM"){print $1, $2-1000, $3+1000}' |sort-bed -) -a <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }' |sort-bed -) | gzip > SNP_in_gencode.vM20.annotation_transcript.bed.gz &


# remove this limit # Keep TRX with more than 10 reads (sum mat/pat F5/F6) (strand specific)
# new requiements: have mat reads >=5 AND pat reads >=5
for Head in HT KD SK
do
  bedtools coverage -s -a <(cat ${Head}_${studyBed}.bed | cut -f 1-6) -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}' > ${Head}_${studyBed}_5+matreads.bed 
  bedtools coverage -s -a ${Head}_${studyBed}_5+matreads.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}'| sort-bed - |uniq > ${Head}_${studyBed}_5+mat_5+patreads_uniq.bed &
done


# the abundance of PolII at each position within the bed file
# sterand specific
# HERE!!!
for Head in HT KD SK
do
  for allele in mat pat
  do
  bedtools coverage -d -s -a ${Head}_${studyBed}_5+mat_5+patreads_uniq.bed -b <(zcat map2ref_1bpbed/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz ) > ${Head}_${allele}_temp.bed  & #|cut -f 8| paste - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > ${Head}_${studyBed}_wSNP_5mat5pat+reads_uniq_${allele}.perBase.bed &
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
