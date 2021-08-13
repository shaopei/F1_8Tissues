cd /workdir/sc2457/F1_Tissues/transcription_level_analysis/TID_TSR_maxTSS

ln -s /workdir/sc2457/F1_Tissues/dREG/Browser/BN_all.dREG.peak.score.bed.gz .
ln -s /workdir/sc2457/F1_Tissues/dREG/Browser/LV_all.dREG.peak.score.bed.gz .

# use maxTSNs
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn/unfiltered_snp.sorted.bed.gz .
for Head in BN LV
do
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn/${Head}_allReads_TSS_maxTSNs.bed .
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn/${Head}_allReads_TSS.bed .
done

# AlleleHMM blocks
ln -s /workdir/sc2457/F1_Tissues/Find_consistent_blocks/HMM_bed . 

Head=BN
for Head in BN LV
do
	# plus and minus strand
zcat ${Head}_all.dREG.peak.score.bed.gz | awk 'BEGIN {OFS="\t"} {print $1,$2,$3, "TID_plus_"NR ,111,"+"} {print $1,$2,$3, "TID_minus_"NR ,111,"-"}' \
> ${Head}_all.dREG.bed


# get allelic read counts for TID, TSS and maxTSN
mkdir toremove
PL=/workdir/sc2457/alleleDB/alleledb_pipeline_mouse
MAPS=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/%s_P.CAST.EiJ_M.C57BL.6J.map
FDR_SIMS=10
FDR_CUTOFF=0.1
BinomialTest(){
  # only use mat-specific or pat-specific reads, ignore IDENTICAL_READ_BED
  # BinomialTest $f $j $MAT_READ_BED $PAT_READ_BED
  # $MAT_READ_BED $PAT_READ_BED need to be sorted!
  # chrX, Y, MT is removed

  f=$1 #bed file for Binomial test eg: BN_allReads_TSS_maxTSNs_SNPs20bp.bed
  j=$2 # name BN_allReads_TSS_maxTSNs_SNPs20bp
  MAT_READ_BED=$3
  PAT_READ_BED=$4

  bedtools coverage -sorted -s -a <(grep -v chrX $f |grep -v chrY | grep -v chrMT |sort-bed -) -b <(zcat ${MAT_READ_BED}| grep -v chrX |grep -v chrY | grep -v chrMT) | cut -f 1-7  > ${j}.mat_cov.bed &
  bedtools coverage -sorted -s -a <(grep -v chrX $f |grep -v chrY | grep -v chrMT |sort-bed -) -b <(zcat ${PAT_READ_BED}| grep -v chrX |grep -v chrY | grep -v chrMT) | cut -f 1-7 > ${j}.pat_cov.bed &
  wait

  # keep every blocks inclding block with at 0 allele-specific read
  paste ${j}.mat_cov.bed ${j}.pat_cov.bed | awk 'BEGIN {OFS="\t"; t="_"} ($2==$9 && $3==$10 && $4==$11 && $6==$13) {print "_", $1, $2, $3, $7,$14, $6}' | \
  awk 'BEGIN{OFS="\t"; t=","} ($5>$6) {print $2, $3, $4, "M"t$5t$6, $5, $6, "_",$7}  ($5<$6) {print $2, $3, $4, "P"t$5t$6, $5, $6, "_",$7}   ($5==$6) {print $2, $3, $4, "S"t$5t$6, $5, $6, "_",$7}  '\
    > ${j}.merged_cov.bed
  wc -l ${j}.mat_cov.bed ${j}.pat_cov.bed  ${j}.merged_cov.bed
  # input format for ${PL}/BinomialTestFor_merged_cov.bed.py: 'chrm','chrmStart', 'chrmEnd', 'hmm_state','mat_allele_count','pat_allele_count','identical_reads_count'...]
  # Col4 mat/pat/Sym states
  # Col5 mat reads count
  # Col6 pat reads count
  #mat  = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[4], skiprows=0)
  #pat = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[5], skiprows=0)

  # output of BinomialTestFor_merged_cov.bed.py:(hmm+BinomialTest) if p-value <= 0.05, remain what it got from hmm (can ne M,P, or S), otherwise S.
  python2 ${PL}/BinomialTestFor_merged_cov.bed.py ${j}.merged_cov.bed ${j}_binomtest.bed
  #R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 0.1 ${j}_binomtest_Rfdr0.1.bed < getCorrectedPValue.R &
  #R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 0.2 ${j}_binomtest_Rfdr0.2.bed < getCorrectedPValue.R &
  #R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 1 ${j}_binomtest_Rfdr1.bed < getCorrectedPValue.R &
  #python2 ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt 
  #awk 'NR==1 { print $0 } NR>1 && ($9+0) <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_binomtest_interestingHets.bed
  mv ${j}.mat_cov.bed ${j}.pat_cov.bed ${j}.merged_cov.bed toremove
}

ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .
bed_dir=map2ref_1bpbed_map5
bowtiebody=bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz
for Head in BN LV #HT  SK  SP  KD  GI  ST
do 
#echo $Head
zcat ${bed_dir}/${Head}_MB6_all_R1.mat.${bowtiebody} ${bed_dir}/${Head}_PB6_all_R1.mat.${bowtiebody} |grep -v chrX |grep -v chrY | grep -v chrMT| sort-bed --max-mem 10G - |gzip >${Head}_mat_temp.gz  &
zcat ${bed_dir}/${Head}_MB6_all_R1.pat.${bowtiebody} ${bed_dir}/${Head}_PB6_all_R1.pat.${bowtiebody} |grep -v chrX |grep -v chrY | grep -v chrMT| sort-bed --max-mem 10G - |gzip >${Head}_pat_temp.gz  &
done
wait

# Perform BinomialTest one strand at a time, using pooled MB6 and PB6 reads 
# mat = B6, pat=CAST
for Head in BN LV # HT  SK  SP  KD  LV  GI  ST
do 
for body in all.dREG allReads_TSS allReads_TSS_maxTSNs
do
    MAT_READ_BED=${Head}_mat_temp.gz
    PAT_READ_BED=${Head}_pat_temp.gz
    BinomialTest ${Head}_${body}.bed ${Head}_${body} ${MAT_READ_BED} ${PAT_READ_BED} &
done
done
wait


ln -s ../pause/P.CAST.EiJ_M.C57BL.6J_*aternal_all.fa* .
Seq-Inr_B6_CAST_TSN(){
j=$1

 #if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 echo "B6_seq" > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt 
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} (NR>1 && $10=="+") {print substr($1,4)p, $2-1, $3, $4,$9,$10} (NR>1 && $10=="-") {print substr($1,4)p, $2, $3+1, $4,$9,$10}')  | grep -v \> >> ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 echo "CAST_seq" > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt 
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} (NR>1 && $10=="+") {print substr($1,4)p, $2-1, $3, $4,$9,$10} (NR>1 && $10=="-") {print substr($1,4)p, $2, $3+1, $4,$9,$10}')   | grep -v \> >> ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 mv ${j}.bed ${j}_temp
 cat ${j}_temp | awk 'BEGIN {OFS="\t"} NR==1 {print $0, "strand"} NR>1 {print $0}'  > ${j}.bed
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_Inr_mat_patSeq.bed 
# cat ${j}_+-${d}_mat_patSeq.bed |awk '{OFS="\t"} NR>1 {print $0}' | awk '{OFS="\t"} (substr($4,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($4,1,1),$11, $12} 
# (substr($4,1,1)=="P") {print  $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($4,1,1), $12, $11}' > ${j}_Inr_High_LowAlleleSeq.bed 
#fi
}

for Head in BN LV #HT  SK  SP  KD  GI  ST
do
  j=${Head}_allReads_TSS_maxTSNs_binomtest
#Seq_High_Low_TSN ${Head} ${k}_binomtest_interestingHets $d
Seq-Inr_B6_CAST_TSN ${j}
done
done

for Head in BN LV #HT  SK  SP  KD  GI  ST
do
# keep TID with at least 5 allelic reads
cat ${Head}_all.dREG_binomtest.bed | awk 'BEGIN {OFS="\t"} (NR==1) {print $1,$2,$3,$4,$5, "strand", $6, $7, $9}( NR>1 && $6+$7 >=5) {print $1,$2,$3,$4,$5, $10, $6, $7, $9}' \
> ${Head}_all.dREG_binomtest_5+reads.bed
# keep TSS with at least 1 allelic reads
cat ${Head}_allReads_TSS_binomtest.bed | awk 'BEGIN {OFS="\t"} (NR==1) {print $1,$2,$3,$4,$5, "strand", $6, $7, $9}( NR>1 && $6+$7 >=1) {print $1,$2,$3,$4,$5, $10, $6, $7, $9}'\
 > ${Head}_allReads_TSS_binomtest_1+read.bed
# TSSs within each TID
bedtools intersect -s -wo -a ${Head}_allReads_TSS_binomtest_1+read.bed -b ${Head}_all.dREG_binomtest_5+reads.bed  \
> ${Head}_TSS_TID.bed

# keep maxTSN with at least 1 allelic reads
cat ${Head}_allReads_TSS_maxTSNs_binomtest_Inr_mat_patSeq.bed | awk 'BEGIN {OFS="\t"} (NR==1) {print $1,$2,$3,$4,$5, "strand", $6, $7, $9, $11, $12}( NR>1 && $6+$7 >=1) {print $1,$2,$3,$4,$5, $10, $6, $7, $9, $11, $12}' \
> ${Head}_allReads_TSS_maxTSNs_binomtest_Inr_mat_patSeq_1+read.bed

# maxTSN within TSS
bedtools intersect -s -wo -a ${Head}_allReads_TSS_maxTSNs_binomtest_Inr_mat_patSeq_1+read.bed -b ${Head}_TSS_TID.bed  \
> ${Head}_maxTSN_TSS_TID.bed
# 1-11 maxTSN, 12-20 TSS, 21-29 TID

# maxTSN overlap with SNPs
bedtools intersect -wo -a ${Head}_maxTSN_TSS_TID.bed -b <(zcat unfiltered_snp.sorted.bed.gz) \
>  ${Head}_maxTSN_TSS_TID_withSNP.bed

# pick TSS(A) with high allele maxTSN CA, low allele nonCA.
cat ${Head}_maxTSN_TSS_TID_withSNP.bed | awk 'BEGIN{OFS="\t"} 
(substr($4,1,1)=="M" && $10 == "CA") {print $0}
(substr($4,1,1)=="P" && $11 == "CA") {print $0}' > ${Head}_maxTSN_TSS_TID_withSNP_temp1.bed

# identify other TSS within TID with TSS(A)
# use the TID(B) contain TSS(A) to identify all TSS(C) within th TID(B)(21-29)
bedtools intersect -s -wo -a ${Head}_allReads_TSS_binomtest_1+read.bed -b <(cat ${Head}_maxTSN_TSS_TID_withSNP_temp1.bed| awk 'BEGIN {OFS="\t"} {print $21,$22, $23, $24, $25, $26, "maxTSN", $4}' | sort-bed - |uniq) \
> ${Head}_maxTSN_TSS_TID_withSNP_temp2.bed
# remove TSS(A) from TSS(C)
bedtools intersect -s -v -a ${Head}_maxTSN_TSS_TID_withSNP_temp2.bed -b <(cat ${Head}_maxTSN_TSS_TID_withSNP_temp1.bed|cut -f 12-20 | sort-bed - |uniq) | uniq \
> ${Head}_maxTSN_TSS_TID_withSNP_temp3.bed
# 1-9 TSS, 10-15 TID, 17 maxTSN with SNP in TSS(A)


# control group : TSS in TID NOT cootain SNPs at any maxTSN
bedtools intersect -v -a ${Head}_maxTSN_TSS_TID.bed -b <(cat ${Head}_maxTSN_TSS_TID_withSNP.bed | cut -f 21-26) \
>  ${Head}_maxTSN_TSS_TID_OutsideTIDwithSNPmaxTSN.bed
# 1-11 maxTSN, 12-20 TSS, 21-29 TID

#HERE

done

TID_TSR_maxTSS_202108.R