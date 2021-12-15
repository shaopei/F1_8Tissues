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

ln -s ../domains_cluster_more_than_chance_or_not_tunit_protein/getFDR.R .
BinomialTest_IDE(){
  # BinomialTest $f $j $MAT_READ_BED $PAT_READ_BED
  # Need to seperate plus and minus strand before using this function. This function ignore strandness in the output
  f=$1  #bed file containing region to perform test
  j=$2  # prefix for output
  MAT_READ_BED=$3
  PAT_READ_BED=$4
  IDE_READ_BED=$5

  bedtools coverage -s -a <(cat $f | LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}') -b <(zcat ${MAT_READ_BED} | grep -v chrX |grep -v chrY | grep -v chrMT) -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3t$6, $1,$2,$3,$7,$8,$9,$10, $6}' > ${j}.mat_cov.bed &
  bedtools coverage -s -a <(cat $f | LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}') -b <(zcat ${PAT_READ_BED} | grep -v chrX |grep -v chrY | grep -v chrMT) -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3t$6, $1,$2,$3,$7,$8,$9,$10, $6}' > ${j}.pat_cov.bed &
  bedtools coverage -s -a <(cat $f | LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}') -b <(zcat ${IDE_READ_BED} | grep -v chrX |grep -v chrY | grep -v chrMT) -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3t$6, $1,$2,$3,$7,$8,$9,$10, $6}' > ${j}.ide_cov.bed &
  wait

  # filter the block and only keep block with at lease 0 allele-specific read
  join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.5,1.9 ${j}.mat_cov.bed ${j}.pat_cov.bed > ${j}.merged_cov.temp
  join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,1.7 ${j}.merged_cov.temp ${j}.ide_cov.bed \
  |awk 'BEGIN{OFS="\t"; t=","} ($5>$6) {print $2, $3, $4, "M"t$5t$6, $5, $6, $7, $8 }  ($5<$6) {print $2, $3, $4, "P"t$5t$6, $5, $6, $7, $8 }  ($5==$6) {print $2, $3, $4, "S"t$5t$6, $5, $6, $7, $8 } '  \
  | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${j}.merged_cov.bed

  mv ${j}.mat_cov.bed ${j}.pat_cov.bed toremove
  # output of BinomialTestFor_merged_cov.bed.py:(hmm+BinomialTest) if p-value <= 0.05, remain what it got from hmm (can ne M,P, or S), otherwise S.
  python2 ${PL}/BinomialTestFor_merged_cov.bed.py ${j}.merged_cov.bed ${j}_binomtest.bed
  mv ${j}.merged_cov.bed toremove
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed ${j}_binomtest_fdr.bed < getFDR.R
  # python2 ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt &
  #  awk 'NR==1 { print $0 } NR>1 && $9 <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_interestingHets.bed
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

ln -s ../domains_cluster_more_than_chance_or_not_tunit_protein/*_READ_BED.temp.gz .


# Perform BinomialTest two strand at a time, using pooled MB6 and PB6 reads 
# mat = B6, pat=CAST
for Head in BN LV # HT  SK  SP  KD  LV  GI  ST
do 
for body in all.dREG allReads_TSS allReads_TSS_maxTSNs
do
    MAT_READ_BED=${Head}_MAT_READ_BED.temp.gz
    PAT_READ_BED=${Head}_PAT_READ_BED.temp.gz
    IDE_READ_BED=${Head}_IDE_READ_BED.temp.gz
    BinomialTest_IDE ${Head}_${body}.bed ${Head}_${body} ${MAT_READ_BED} ${PAT_READ_BED} ${IDE_READ_BED} &
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

# maxTSN overlap with SNPs, # TSSs within each TID
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


# control group : TSS in TID NOT contain SNPs at any maxTSN
bedtools intersect -v -a ${Head}_maxTSN_TSS_TID.bed -b <(cat ${Head}_maxTSN_TSS_TID_withSNP.bed | cut -f 21-26) \
>  ${Head}_maxTSN_TSS_TID_OutsideTIDwithSNPmaxTSN.bed
# 1-11 maxTSN, 12-20 TSS, 21-29 TID

#HERE

done

TID_TSR_maxTSS_202108.R

### identify sequence upstream of maxTSN -35 to -20 ,(TATA-box -31 to -24)
Seq-TATA_B6_CAST_TSN(){
j=$1

 #if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 echo "B6_seq" > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt 
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} (NR>1 && $10=="+") {print substr($1,4)p, $2-35, $3-20, $4,$9,$10} (NR>1 && $10=="-") {print substr($1,4)p, $2+20, $3+35, $4,$9,$10}')  | grep -v \> >> ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 echo "CAST_seq" > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt 
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} (NR>1 && $10=="+") {print substr($1,4)p, $2-35, $3-20, $4,$9,$10} (NR>1 && $10=="-") {print substr($1,4)p, $2+20, $3+35, $4,$9,$10}')   | grep -v \> >> ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_TATA_mat_patSeq.bed 
# cat ${j}_+-${d}_mat_patSeq.bed |awk '{OFS="\t"} NR>1 {print $0}' | awk '{OFS="\t"} (substr($4,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($4,1,1),$11, $12} 
# (substr($4,1,1)=="P") {print  $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($4,1,1), $12, $11}' > ${j}_Inr_High_LowAlleleSeq.bed 
#fi
}

for Head in BN LV #HT  SK  SP  KD  GI  ST
do
  j=${Head}_allReads_TSS_maxTSNs_binomtest
  # the sequence upstream of maxTSN
  Seq-TATA_B6_CAST_TSN ${j}
done

# maxTSN within TSS 1+ reads, TID 5+ reads
bedtools intersect -s -wo -a <(cat ${Head}_allReads_TSS_maxTSNs_binomtest_TATA_mat_patSeq.bed| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5, $10,$6,$7,$8, $9, $11, $12}') -b ${Head}_TSS_TID.bed  \
> ${Head}_maxTSN.TATA_TSS_TID.bed
#chrm   chrmStart       chrmEnd hmm_state       hmm+BinomialTest   strand      mat_allele_count        pat_allele_count        identical_reads_count   Binom_p_value        B6_seq  CAST_seq
# 1-12 maxTSN, 13-21 TSS, 22-30 TID

# maxTSN with TSS with different AT content in maxTSN upstream TATA region
# ${Head}_maxTSN.TATA_TSS_TID_diffAT.bed from TID_TSR_maxTSS_202108.R
# identify other TSS within TID with TSS(A)
# use the TID(B) contain TSS(A) to identify all TSS(C) within th TID(B)(22-27)
bedtools intersect -s -wo -a ${Head}_allReads_TSS_binomtest_1+read.bed -b <(cat ${Head}_maxTSN.TATA_TSS_TID_diffAT.bed| awk 'BEGIN {OFS="\t"} {print $22, $23, $24, $25, $26, $27, "maxTSN", $4, $11, $12}' | sort-bed - |uniq) \
> ${Head}_maxTSN.TATA_TSS_TID_diffAT_TSSC.bed
# remove TSS(A) from TSS(C)
bedtools intersect -s -v -a ${Head}_maxTSN.TATA_TSS_TID_diffAT_TSSC.bed -b <(cat ${Head}_maxTSN.TATA_TSS_TID_diffAT.bed|cut -f 13-18 | sort-bed - |uniq) | uniq \
> ${Head}_maxTSN.TATA_TSS_TID_diffAT_TSSC-A.bed
# 1-9 TSS, 10-15 TID, 17 maxTSN with different AT content in TSS(A)


### overlap of TSCs with (A)allelic difference in abundance; (B)allelic difference in shape, and (C)inside AlleleHMM blocks


for Head in BN LV #HT  SK  SP  KD  GI  ST
do
  # (A)allelic difference in abundance
  # j=${Head}_allReads_TSS_binomtest_1+read   
  R --vanilla --slave --args $(pwd) ${Head}_allReads_TSS_binomtest_1+read.bed ${Head}_allReads_TSS_binomtest_1+read_fdr.bed < getFDR.R
  cat ${Head}_allReads_TSS_binomtest_1+read_fdr.bed | awk 'BEGIN{OFS="\t"} $10+0<=0.1 {print $0}' > ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed

 # (B)allelic difference in shape from F1_TSS_KStest_MultiBaseRunOn_Formanuscript.sh
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSS_KStest_MultiBaseRunOn/${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed
#ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSS_KStest_MultiBaseRunOn/${Head}_allReads_TSS_5mat5pat_uniq_pValue.bed  # this one contain more info indlucing p-value and fdr

# (C)TSCs within AlleleHMM blocks
intersectBed -s -u -a ${Head}_allReads_TSS.bed -b <(cat HMM_bed/${Head}_*HMM*bed) > ${Head}_allReads_TSS_inAlleleHMMBlocks.bed

#intersect of A & B
wc -l ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed
wc -l ${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed
wc -l ${Head}_allReads_TSS_inAlleleHMMBlocks.bed
cat ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed |awk 'BEGIN {OFS="\t";l="_"}{print $1l$2l$3$6}' > ${Head}_temp_A
cat ${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed |awk 'BEGIN {OFS="\t";l="_"}{print $1l$2l$3$6}' > ${Head}_temp_B
cat ${Head}_allReads_TSS_inAlleleHMMBlocks.bed |awk 'BEGIN {OFS="\t";l="_"}{print $1l$2l$3$6}' > ${Head}_temp_C

# A,B
intersectBed -s -wo -a ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed -b ${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed | wc -l 
# A,C
intersectBed -s -wo -a ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed -b ${Head}_allReads_TSS_inAlleleHMMBlocks.bed | wc -l 
# B,C
intersectBed -s -wo -a ${Head}_allReads_TSS_inAlleleHMMBlocks.bed -b ${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed | wc -l 
#A.B.C
intersectBed -s -wo -a ${Head}_allReads_TSS_inAlleleHMMBlocks.bed -b <(intersectBed -s -u -a ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed -b ${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed) |wc -l

done

### plot venn diagram in R
library(eulerr)
set.seed(1)
#BN
combo <- c(A = 3225, B = 1093, C = 2324, "A&B" = 476, "A&C" = 410, "B&C" = 126, "A&B&C"=94)
plot(euler(combo))
#LV
combo <- c(A = 7227, B = 1482, C = 10113, "A&B" = 745, "A&C" = 2133, "B&C" = 313, "A&B&C"=242)
plot(euler(combo))

###
#R CMD INSTALL /Users/sc2457/Downloads/BioVenn_1.1.3.tar.gz
library(BioVenn)
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/TID_TSR_maxTSS")
Head="BN"
a=read.table(paste(Head,"_temp_A", sep=""))
b=read.table(paste(Head,"_temp_B", sep=""))
c=read.table(paste(Head,"_temp_C", sep=""))
biovenn <- draw.venn(a$V1, b$V1, c$V1, subtitle=Head, nrtype="abs")
biovenn <- draw.venn(a$V1, b$V1, c$V1, subtitle="Example diagram 1", nrtype="abs")


### Changes in the abundance of Pol II in TSCs were enriched in changes in AlleleHMM blocks (give results)
# No allelic difference in abundance (fdr0.9) inside alleleHMM blocks
# continue...
# A,C
Head=BN
wc -l ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed
intersectBed -s -wo -a ${Head}_allReads_TSS_binomtest_1+read_fdr0.1.bed -b ${Head}_allReads_TSS_inAlleleHMMBlocks.bed | wc -l 
cat ${Head}_allReads_TSS_binomtest_1+read_fdr.bed | awk 'BEGIN {OFS="\t"} ($10+0>0.9){print $0}' > ${Head}_allReads_TSS_binomtest_1+read_fdr0.9.bed
wc -l ${Head}_allReads_TSS_binomtest_1+read_fdr0.9.bed
intersectBed -s -wo -a ${Head}_allReads_TSS_binomtest_1+read_fdr0.9.bed -b ${Head}_allReads_TSS_inAlleleHMMBlocks.bed | wc -l 
# in R
fisher.test(matrix(c(dim(df_withoutAT)[1],sum(df_withoutAT$pair_fdr <= 0.1), 
                     dim(unique(df_withAT[,1:6]))[1],dim(unique(df_withAT[df_withAT$pair_fdr <= 0.1, 1:6]))[1])
             , 2,2))

####
# identify TBP binding motif using rtdbsdb
# Merge M00216_2.00 and M09433_2.00
# Use the maxTSN to identify the TATA box upstream -20~-35bp

cat mm10-TBP_M00216_2.00.bed | awk 'BEGIN{OFS="\t"}(NR>1){print $0}') > mm10-TBP_M00216andM09433.bed
cat mm10-TBP_M09433_2.00.bed | awk 'BEGIN{OFS="\t"}(NR>1){print $0}') >> mm10-TBP_M00216andM09433.bed
sort-bed mm10-TBP_M00216andM09433.bed > mm10-TBP_M00216andM09433_sorted.bed
# Use the maxTSN to identify the TATA box upstream -20~-35bp
Head=LV  #BN
bedtools closest -s -D a -id -a <(cat ${Head}_allReads_TSS_maxTSNs_binomtest.bed | awk 'BEGIN{OFS="\t"}(NR>1){print $1,$2,$3,$4,$5, $10,$6,$7,$8, $9}' | sort-bed - ) \
-b mm10-TBP_M00216andM09433_sorted.bed  | awk 'BEGIN{OFS="\t"}($19 >= -35 && $19 <=-20){print $0}'> ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATA.bed
# col1-10 maxTSN; col 11-18 tatabox
#   608 BN_allReads_TSS_maxTSNs_binomtest_closestTATA.bed
#   297 LV_allReads_TSS_maxTSNs_binomtest_closestTATA.bed

# seperate to with SNP in tatabox, and without SNP in tatabox
# identify maxTSN with SNPs in the TATAbox -35~-20 upstream                                      # col1-10 maxTSN; col 11-18 tatabox
bedtools intersect -sorted -wo -a <(cat ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATA.bed |awk 'BEGIN{OFS="\t"}{print $11,$12,$13,$14,$15,$16,$17,$18,$1,$2,$3,$4,$5,$6,$7,$8,$9, $10, $19}')\
 -b <(zcat unfiltered_snp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"}{print $9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$1,$2,$3,$4,$5,$6,$7,$8, $19, $20,$21,$22,$23}' \
|sort-bed - | uniq > ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATAwithSNP.bed
 # ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATAwithSNP.bed: col1-10 maxTSN; col 11-18 tatabox; col19 dis(maxTSN, TATA); col20-23 SNPs info
#  28 BN_allReads_TSS_maxTSNs_binomtest_closestTATAwithSNP.bed
#  20 LV_allReads_TSS_maxTSNs_binomtest_closestTATAwithSNP.bed


# identify maxTSN withOUT SNPs in the TATAbox -35~-20 upstream 
bedtools intersect -sorted -v -a <(cat ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATA.bed |awk 'BEGIN{OFS="\t"}{print $11,$12,$13,$14,$15,$16,$17,$18,$1,$2,$3,$4,$5,$6,$7,$8,$9, $10, $19}')\
 -b <(zcat unfiltered_snp.sorted.bed.gz) | awk 'BEGIN{OFS="\t"}{print $9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$1,$2,$3,$4,$5,$6,$7,$8, $19}' \
|sort-bed - | uniq > ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATAwithoutSNP.bed
 # ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATAwithoutSNP.bed: col1-10 maxTSN; col 11-18 tatabox; col19 dis(maxTSN, TATA); 

# identify the TSS(A) and TID(B) contain maxTSN with SNPs in TATAbox
bedtools intersect -s -wo -a <(cat ${Head}_allReads_TSS_maxTSNs_binomtest_closestTATAwithSNP.bed |cut -f 1-16) -b ${Head}_TSS_TID.bed  \
> ${Head}_maxTSN.closestTATAwithSNP_TSS_TID.bed
#col1-10 maxTSN, 11-16 tatabox, 17-25 TSS(A), 26-31 TID(B)


# identify other TSS within TID(B) with TSS(A)
# use the TID(B) contain TSS(A) to identify all TSS(C) within th TID(B)(20-25)
bedtools intersect -s -wo -a ${Head}_allReads_TSS_binomtest_1+read.bed \
-b <(cat ${Head}_maxTSN.closestTATAwithSNP_TSS_TID.bed| awk 'BEGIN {OFS="\t"} {print $26,$27, $28, $29, $30, $31,  "maxTSN", $4, $7, $8, $9, "tatabox", $11,$12,$13,$14, $15,$16}' | sort-bed - |uniq) \
> ${Head}_maxTSN.closestTATAwithSNP_TSS_TID_TSSC.bed
# remove TSS(A) from TSS(C)
bedtools intersect -s -v -a ${Head}_maxTSN.closestTATAwithSNP_TSS_TID_TSSC.bed -b <(cat ${Head}_maxTSN.closestTATAwithSNP_TSS_TID.bed|cut -f 17-22 | sort-bed - |uniq) | uniq \
> ${Head}_maxTSN.closestTATAwithSNP_TSS_TID_TSSC-A.bed
# 1-9 TSS, 10-15 TID, 16-20 maxTSN with SNPs in tatabox in TSS(A), 21-27 tatabox


# get the tatabox sequence from ${Head}_maxTSN.closestTATAwithSNP_TSS_TID_TSSC-A.bed
j=${Head}_maxTSN.closestTATAwithSNP_TSS_TID_TSSC-A
 #if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed | cut -f 22-27 |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2, $3, $4,$5,$6 }')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed | cut -f 22-27 | awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2, $3, $4,$5,$6 }')| grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_TATA_mat_patSeq.bed 




### 
# identify pontential Tata box region -35~-20 upstream of maxTSN
# give fasta files and bed regions to rtdbsdb, 
# use tfbs.scanTFsite to get the score of TBP binding motif in the given bed regions
# compare the score between B6 and CAST

for Head in BN LV #HT  SK  SP  KD  GI  ST
do
  j=${Head}_allReads_TSS_maxTSNs_binomtest
 cat ${j}.bed |awk '{OFS="\t";p="_maternal"} (NR>1 && $10=="+") {print substr($1,4)p, $2-35, $3-20, $4,$9,$10} (NR>1 && $10=="-") {print substr($1,4)p, $2+20, $3+35, $4,$9,$10}' > ${j}_-35to-20_mat.bed
 cat ${j}.bed |awk '{OFS="\t";p="_paternal"} (NR>1 && $10=="+") {print substr($1,4)p, $2-35, $3-20, $4,$9,$10} (NR>1 && $10=="-") {print substr($1,4)p, $2+20, $3+35, $4,$9,$10}' > ${j}_-35to-20_pat.bed

done

# rtdbsdb_mm10_TATA.R

# Use the upstream -20~-35bp of maxTSN to represent the pontentail TATA box region 
# -a -35~-20 upstream of maxTSN
# -b motif binding sites from rtdbsdb_mm10_TATA.R, score of all regions with binding motif
Head=LV
 cat ${Head}_allReads_TSS_maxTSNs_binomtest.bed |awk '{OFS="\t"} (NR>1 && $10=="+") {print $1, $2-35, $3-20, $4,$9,$10, $6,$7,$8, $1, $2, $3} (NR>1 && $10=="-") {print $1, $2+20, $3+35, $4,$9,$10, $6,$7,$8, $1, $2, $3}' > test1.bed
 # chrm   chrmStart  chrmEnd hmm_state   Binom_p_value   strand   mat_allele_count   pat_allele_count   identical_reads_count
 bedtools intersect -wo -a test1.bed \
 -b <(cat ${Head}_allReads_TSS_maxTSNs_binomtest_motifM00216.bed |awk '{OFS="\t"} (NR>1) {print "chr"$0}') \
 > ${Head}_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM00216.bed

 bedtools intersect -wo -a test1.bed \
 -b <(cat ${Head}_allReads_TSS_maxTSNs_binomtest_motifM09433.bed |awk '{OFS="\t"} (NR>1) {print "chr"$0}') \
 > ${Head}_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM09433.bed



# -a chrm(TATA region)   chrmStart  chrmEnd hmm_state   Binom_p_value   strand   mat_allele_count   pat_allele_count   identical_reads_count chrm(maxTSN)   chrmStart  chrmEnd
# -b chrmo chromStart  chromEnd        score.b6        score.cast      strand  score.ave

# get the maxScore of avrage(score.B6, score.CAST) from the intersect motif bind sites # only use the first, if tie
python2 getMaxScore_TFbindingmotif.py ${Head}_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM00216.bed ${Head}_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM00216_maxScore.bed
python2 getMaxScore_TFbindingmotif.py ${Head}_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM09433.bed ${Head}_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM09433_maxScore.bed





