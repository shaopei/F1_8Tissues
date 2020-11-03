cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSN_ShootingGallery/WithAlleleHMMFilter_strandBigCorrected

# use maxTSNs
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn/*_allReads_TSS_maxTSNs.bed .
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn/unfiltered_snp.sorted.bed.gz .
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn/*_allReads_TSS.bed .


#Head=BN
#Head=LV

# wc -l BN_allReads_TSS.bed
#78652 BN_allReads_TSS.bed

# remove the TSS overlap with AlleleHMM blocks
# remove chrX
ln -s /workdir/sc2457/F1_Tissues/Find_consistent_blocks/HMM_bed . # AlleleHMM blocks

for Head in BN LV
do
bedtools intersect -v -a <(grep -v chrX ${Head}_allReads_TSS.bed) -b <(cat HMM_bed/combined_cross/${Head}_HMM_*.bed) > ${Head}_allReads_TSS_NotInAlleleHMMBlocks.bed

# maxTSN with SNPs at initiation motif
bedtools closest -D a -id -a <(sort-bed ${Head}_allReads_TSS_maxTSNs.bed) -b <(zcat unfiltered_snp.sorted.bed.gz) | awk 'BEGIN {OFS="\t"; t="_"} ($11==-1 || $11==0) {print $1, $2, $3, $4, $5, $6}'  > ${Head}_allReads_TSS_maxTSNs_SNP.bed
# -id	Ignore features in B that are downstream of features

# maxTSN with SNPs at initiation motif within the TSS NOT overlap with AlleleHMM blocks
bedtools intersect -s -a ${Head}_allReads_TSS_maxTSNs_SNP.bed -b ${Head}_allReads_TSS_NotInAlleleHMMBlocks.bed > ${Head}_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks.bed
#body=allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks


# maxTSN with or without SNPs at initiation motif within the TSS NOT overlap with AlleleHMM blocks
bedtools intersect -s -a ${Head}_allReads_TSS_maxTSNs.bed -b ${Head}_allReads_TSS_NotInAlleleHMMBlocks.bed > ${Head}_allReads_TSS_maxTSNs_TSSNotInAlleleHMMBlocks.bed

# maxTSN WITHOUT SNPs at initiation motif within the TSS NOT overlap with AlleleHMM blocks
bedtools intersect -s -a <(bedtools intersect -v -s -a ${Head}_allReads_TSS_maxTSNs.bed -b ${Head}_allReads_TSS_maxTSNs_SNP.bed) -b ${Head}_allReads_TSS_NotInAlleleHMMBlocks.bed > ${Head}_allReads_TSS_maxTSNs_NoSNP_TSSNotInAlleleHMMBlocks.bed


#body=allReads_TSS_maxTSNs_TSSNotInAlleleHMMBlocks
 wc -l ${Head}_allReads_TSS_maxTSNs*_TSSNotInAlleleHMMBlocks.bed
   # 83088 BN_allReads_TSS_maxTSNs.bed
   #   612 BN_allReads_TSS_maxTSNs_SNP.bed
   #   520 BN_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks.bed
   # 76398 BN_allReads_TSS_maxTSNs_TSSNotInAlleleHMMBlocks.bed
   # 97577 LV_allReads_TSS_maxTSNs.bed
   #   935 LV_allReads_TSS_maxTSNs_SNP.bed
   #   673 LV_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks.bed
   # 77574 LV_allReads_TSS_maxTSNs_TSSNotInAlleleHMMBlocks.bed

  # 75878 BN_allReads_TSS_maxTSNs_NoSNP_TSSNotInAlleleHMMBlocks.bed
  #   520 BN_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks.bed
  # 76398 BN_allReads_TSS_maxTSNs_TSSNotInAlleleHMMBlocks.bed
  
  # 76901 LV_allReads_TSS_maxTSNs_NoSNP_TSSNotInAlleleHMMBlocks.bed
  #   673 LV_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks.bed
  # 77574 LV_allReads_TSS_maxTSNs_TSSNotInAlleleHMMBlocks.bed


done

mkdir toremove
PL=/workdir/sc2457/alleleDB/alleledb_pipeline_mouse
MAPS=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/%s_P.CAST.EiJ_M.C57BL.6J.map
FDR_SIMS=10
FDR_CUTOFF=0.1

BinomialTest_TSN(){
  # only use mat-specific or pat-specific reads, ignore IDENTICAL_READ_BED
  # take stradness into account but only process one strand (plus or minus) at a time
  # BinomialTest $f $j $MAT_READ_BED $PAT_READ_BED
  # $MAT_READ_BED $PAT_READ_BED need to be sorted!
  # chrX is removed

  f=$1 #bed file for Binomial test eg: BN_allReads_TSS_maxTSNs_SNPs20bp.bed
  j=$2 # name BN_allReads_TSS_maxTSNs_SNPs20bp
  MAT_READ_BED=$3
  PAT_READ_BED=$4

  bedtools coverage -sorted -s -a <(grep -v chrX $f |sort-bed -) -b <(zcat ${MAT_READ_BED}) | cut -f 1-7  > ${j}.mat_cov.bed &
  bedtools coverage -sorted -s -a <(grep -v chrX $f |sort-bed -) -b <(zcat ${PAT_READ_BED}) | cut -f 1-7 > ${j}.pat_cov.bed &
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
  python2 ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt 
  awk 'NR==1 { print $0 } NR>1 && ($9+0) <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_binomtest_interestingHets.bed
  mv ${j}.mat_cov.bed ${j}.pat_cov.bed ${j}.merged_cov.bed toremove
}


# Pool MB6 and PB6 reads together. mat is B6 reads, pat is Cast reads
# remove chrX
ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .
bed_dir=map2ref_1bpbed_map5
bowtiebody=bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz
for Head in BN LV #HT  SK  SP  KD  GI  ST
do 
#echo $Head
zcat ${bed_dir}/${Head}_MB6_all_R1.mat.${bowtiebody} ${bed_dir}/${Head}_PB6_all_R1.mat.${bowtiebody} |grep -v chrX | sort-bed --max-mem 10G - |gzip >${Head}_mat_temp.gz  &
zcat ${bed_dir}/${Head}_MB6_all_R1.pat.${bowtiebody} ${bed_dir}/${Head}_PB6_all_R1.pat.${bowtiebody} |grep -v chrX | sort-bed --max-mem 10G - |gzip >${Head}_pat_temp.gz  &
done
wait

# Perform BinomialTest one strand at a time, using pooled MB6 and PB6 reads 
for Head in BN LV # HT  SK  SP  KD  LV  GI  ST
do 
for body in allReads_TSS_maxTSNs allReads_TSS_maxTSNs_SNP allReads_TSS_maxTSNs_NoSNP_TSSNotInAlleleHMMBlocks allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks
do
    MAT_READ_BED=${Head}_mat_temp.gz
    PAT_READ_BED=${Head}_pat_temp.gz
    BinomialTest_TSN ${Head}_${body}.bed ${Head}_${body} ${MAT_READ_BED} ${PAT_READ_BED} &
done
done
wait

# get sequence for High and Low Allele (defind by TSN expression level)
ln -s ../P.CAST.EiJ_M.C57BL.6J_*aternal_all.fa* .
Seq_High_Low_TSN(){
 Head=$1
 j=$2
 d=$3
 #step=$4

 #if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 echo "mat_seq" > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt 
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} NR>1 {print substr($1,4)p, $2-d, $3+d, $4,$9,$10}')  | grep -v \> >> ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 echo "pat_seq" > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt 
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} NR>1 {print substr($1,4)p, $2-d, $3+d, $4,$9,$10}')  | grep -v \> >> ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_mat_patSeq.bed 
 cat ${j}_+-${d}_mat_patSeq.bed |awk '{OFS="\t"} NR>1 {print $0}' | awk '{OFS="\t"} (substr($4,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($4,1,1),$11, $12} 
 (substr($4,1,1)=="P") {print  $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($4,1,1), $12, $11}' > ${j}_+-${d}_High_LowAlleleSeq.bed 
#fi

#R --vanilla --slave --args ${j}_+-${d}_High_LowAlleleSeq.bed  ${Head} ${step} ${d} < getGC_content_HighLowAllele.R 
}

for Head in BN LV #HT  SK  SP  KD  GI  ST
do
for body in allReads_TSS_maxTSNs allReads_TSS_maxTSNs_SNP  allReads_TSS_maxTSNs_NoSNP_TSSNotInAlleleHMMBlocks allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks
do
	d=50
  k=${Head}_${body}
#Seq_High_Low_TSN ${Head} ${k}_binomtest_interestingHets $d
Seq_High_Low_TSN ${Head} ${k}_binomtest $d
done
done

# seperate into TSS in Sinlge Base driven and in multiple base driven
for body in allReads_TSS_maxTSNs_NoSNP_TSSNotInAlleleHMMBlocks allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks
do
for Head in BN LV
do
for asTSS in SingleBaseDriven MultipleBaseDriven
do
    ln -s ../${Head}_allReads_TSS_5mat5pat_uniq_AsTSS${asTSS}.bed .
bedtools intersect -u -a <(cat ${Head}_${body}_binomtest_+-50_High_LowAlleleSeq.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5, $10, $0}') \
-b ${Head}_allReads_TSS_5mat5pat_uniq_AsTSS${asTSS}.bed |cut -f 7- > ${Head}_${body}_binomtest_+-50_High_LowAlleleSeq_AsTSS${asTSS}.bed
done
done
done


R script TSN_HighLowAllele_forManuscript.R