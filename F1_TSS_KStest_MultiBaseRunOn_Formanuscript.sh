cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSS_KStest_MultiBaseRunOn

# TSS identify see https://github.com/shaopei/F1_8Tissues/blob/master/F1_TSN_identifyTSS_MultiBaseRunOn.sh
# TSS are strand-specific
ln -s ../identifyTSS_MultiBaseRunOn/*_allReads_TSS.bed .
studyBed=allReads_TSS

# Use TSS sites with  mat reads >=5 AND pat reads >=5 (strand specific) 
# remove TSS in chrX
ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -s -a <(bedtools coverage -sorted -a <(grep -v chrX ${Head}_${studyBed}.bed) -b <(zcat map2ref_1bpbed_map5/${Head}_*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz | sort-bed - --max-mem 10G ) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5,$6}') \
   -b <(zcat map2ref_1bpbed_map5/${Head}_*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz | sort-bed - --max-mem 10G ) | awk 'BEGIN{OFS="\t"} ($7 >=5){print $1,$2,$3,$4,$5, $6}'\
   | sort-bed - --max-mem 10G |uniq > ${Head}_${studyBed}_5mat5pat_uniq.bed &
done
wait

[sc2457@cbsudanko TSS_KStest_MultiBaseRunOn]$ wc -l *bed
   78652 BN_allReads_TSS.bed
   10991 BN_allReads_TSS_5mat5pat_uniq.bed
   41125 GI_allReads_TSS.bed
    7180 GI_allReads_TSS_5mat5pat_uniq.bed
   42521 HT_allReads_TSS.bed
    6747 HT_allReads_TSS_5mat5pat_uniq.bed
   39004 KD_allReads_TSS.bed
    6662 KD_allReads_TSS_5mat5pat_uniq.bed
   91606 LV_allReads_TSS.bed
   17595 LV_allReads_TSS_5mat5pat_uniq.bed
   32363 SK_allReads_TSS.bed
    4568 SK_allReads_TSS_5mat5pat_uniq.bed
   53686 SP_allReads_TSS.bed
    9280 SP_allReads_TSS_5mat5pat_uniq.bed
   34172 ST_allReads_TSS.bed
    5223 ST_allReads_TSS_5mat5pat_uniq.bed
  481375 total



# identify the abundance of PolII at each position
# strand specific
# use all reads (not just allelic reads)
wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz |sort-bed - --max-mem 10G  ) > ${Head}_allReads_temp0.bed &
done

# identify the abundance of PolII at each position
# strand specific
# use ONLY allelic reads
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
  bedtools coverage -sorted -d -s -a ${Head}_${studyBed}_5mat5pat_uniq.bed -b <(zcat map2ref_1bpbed_map5/${Head}_*_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed - --max-mem 10G) > ${Head}_${allele}_temp0.bed &
done
done
wait 

# paste the all reads  and alleleic reads count into a table
# only keep base with at least b all reads ($8 >= b)
b=5
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
    paste ${Head}_allReads_temp0.bed ${Head}_${allele}_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"} ($2==$10 && $3==$11 && $8+0 >= b ) {print $0}' |cut -f 9- > ${Head}_${allele}_temp.bed &
  done
done

wait
ln -s ../Generate_vector_input_for_KStest_NodupPlusMinus.py .
ln -s ../KStest_flexible_length.R .
# use python script to generaye input for KS test
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
python2 Generate_vector_input_for_KStest_NodupPlusMinus.py ${Head}_${allele}_temp.bed ${Head}_${studyBed}_5mat5pat_uniq_${allele}.perBase.bed &
done
done
wait

# get p-value for KS test in R
for Head in BN HT  SK  SP  KD  LV  GI  ST 
do
R --vanilla --slave --args $(pwd) ${Head} ${studyBed}_5mat5pat_uniq < KStest_flexible_length.R &
# output ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed bed6, col7 p.value, col8 fdr
done
wait
# exmaine TSS in IGV
f=0.1
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= f){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_fdr${f}.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= f && $6 =="+"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_fdr${f}_plus.bed
  cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"; c=":"; d="-"} ($8+0 <= f && $6 =="-"){print $1, $2, $3, $4c$6c$8,$5,$6 }' |sort-bed - > ${Head}_${studyBed}_5mat5pat_uniq_fdr${f}_minus.bed
done

###
### seperate TSS with some kind of allelic difference into two group: 1. change driven by 1 single base 2. Change driven by more than 1 base
# examine all TSS
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSS_KStest_MultiBaseRunOn
## mask the TSN with biggest difference between the two alleles
# identify the abundance of PolII at each position
# strand specific
# use all reads (not just allelic reads)



# mask the TSN with biggest difference between the two alleles
# for TSS with only 1bp length, it will be removed after the masking
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
paste ${Head}_mat_temp.bed ${Head}_pat_temp.bed  | awk -v b=$b 'BEGIN{OFS="\t"} ($2==$10 && $3==$11 && $7 == $15 ) {print $0}'> ${Head}_mat_pat_temp.bed 
python2 MaskTSNwithMaxAllelicDifference.py ${Head}_mat_pat_temp.bed ${Head}_mat_pat_temp_masked.bed 
cat ${Head}_mat_pat_temp_masked.bed | cut -f 1-8 > ${Head}_mat_temp_masked.bed
cat ${Head}_mat_pat_temp_masked.bed | cut -f 9-16 > ${Head}_pat_temp_masked.bed
done

wait
ln -s ../Generate_vector_input_for_KStest_NodupPlusMinus.py .
ln -s ../KStest_flexible_length.R .
# use python script to generaye input for KS test
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  for allele in mat pat
  do
python2 Generate_vector_input_for_KStest_NodupPlusMinus.py ${Head}_${allele}_temp_masked.bed ${Head}_${studyBed}_5mat5pat_uniq_masked_${allele}.perBase.bed &
done
done
wait

# get p-value for KS test in R
# the input bed region that do not have at least 5 pat AND 5 mat reads will be removed in KStest_flexible_length.R
for Head in BN HT  SK  SP  KD  LV  GI  ST 
do
cat ${Head}_${studyBed}_5mat5pat_uniq_masked_mat.perBase.bed | cut -f 1-6 > ${Head}_${studyBed}_5mat5pat_uniq_masked.bed 
R --vanilla --slave --args $(pwd) ${Head} ${studyBed}_5mat5pat_uniq_masked < KStest_flexible_length.R &
# output ${Head}_${studyBed}_5mat5pat_uniq_masked_pValue.bed bed6, col7 p.value, col8 fdr
done

# make a table to comapre masked and unmasked p-value and fdr
# -a masked -b umasked
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
bedtools intersect -wb -s -a ${Head}_${studyBed}_5mat5pat_uniq_masked_pValue.bed \
-b ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed | awk -v f=$f 'BEGIN{OFS="\t"} {print $9,$10,$11,$12,$13,$14,$7,$8,$15,$16}' \
> ${Head}_${studyBed}_5mat5pat_uniq_maskedVSunmasked_pValue.bed  &
# output is bed6, col7 maksed pavlue,col 8 masked fdr,col 9 unmaksed pvalue, col 10 unmaked fdr
done

# label the location of maxTSS
 ln -s ../identifyTSS_maxTSNs_MultiBaseRunOn/*_allReads_TSS_maxTSNs.bed .
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
 bedtools intersect -wao -s -a ${Head}_${studyBed}_5mat5pat_uniq_maskedVSunmasked_pValue.bed -b ${Head}_allReads_TSS_maxTSNs.bed > ${Head}_${studyBed}_5mat5pat_uniq_maskedVSunmasked_pValue_maxTSNs.bed &
# distribuiton of SNPs in TSS as comparation
 bedtools intersect -wao -s -a <(cat ${Head}_${studyBed}_5mat5pat_uniq_pValue.bed |awk '{OFS="\t"} {print $1, $2,$3,$4, $5,$6, "NA", "NA", $7,$8 }' ) -b ${Head}_allReads_TSS_maxTSNs.bed > ${Head}_${studyBed}_5mat5pat_uniq_pValue_maxTSNs.bed &
done
ln -s ../P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw .
#Rscript getSNPsAbundance.R
Rscript SNPsAbundance_manuscript_figure.R

# Organ=c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")
# s_count <- NULL
# m_count <- NULL
# for(t in Organ){
#   temp = SNPsAbundanceAroundMaxTSNInTSS(maxTSN, d=202, step=4,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE)
#   s_count = c(s_count, temp[[4]][1])
#   m_count = c(m_count, temp [[5]][1])
#   }
# df = data.frame(Organ, s_count, m_count)




#### Determined High/Low allele based on the transcrition level at TSS
### Binomial test of  TSS 
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/TSS_BinomialTest_MultiBaseRunOn
ln -s ../identifyTSS_MultiBaseRunOn/*_allReads_TSS.bed .
ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .
ln -s ../identifyTSS_MultiBaseRunOn/getCorrectedPValue.R .

PL=/workdir/sc2457/alleleDB/alleledb_pipeline_mouse
MAPS=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/%s_P.CAST.EiJ_M.C57BL.6J.map
FDR_SIMS=10
FDR_CUTOFF=0.1

wait_a_second() {
  joblist=($(jobs -p))
    while (( ${#joblist[*]} >= 10 ))
      do
      sleep 1
      joblist=($(jobs -p))
  done
}

BinomialTest_TSS(){
  # only use mat-specific or pat-specific reads, ignore IDENTICAL_READ_BED
  # take stradness into account but only process one strand (plus or minus) at a time
  # BinomialTest $f $j $MAT_READ_BED $PAT_READ_BED
  # $MAT_READ_BED $PAT_READ_BED need to be sorted!
  # chrX is removed

  f=$1 #bed file for Binomial test eg: BN_allReads_TSS.bed
  j=$2 # name BN_allReads_TSS
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
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 0.1 ${j}_binomtest_Rfdr0.1.bed < getCorrectedPValue.R &
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 1 ${j}_binomtest_Rfdr1.bed < getCorrectedPValue.R &
  #python2 ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt 
  #awk 'NR==1 { print $0 } NR>1 && ($9+0) <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_binomtest_interestingHets.bed
  mv ${j}.mat_cov.bed ${j}.pat_cov.bed ${j}.merged_cov.bed toremove
}


# Pool MB6 and PB6 reads together. mat is B6 reads, pat is Cast reads
# remove chrX
bed_dir=map2ref_1bpbed_map5
body=bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz

for Head in BN HT  SK  SP  KD  LV  GI  ST
do 
#echo $Head
zcat ${bed_dir}/${Head}_MB6_all_R1.mat.${body} ${bed_dir}/${Head}_PB6_all_R1.mat.${body} |grep -v chrX | sort-bed --max-mem 10G - |gzip >${Head}_mat_temp.gz  &
zcat ${bed_dir}/${Head}_MB6_all_R1.pat.${body} ${bed_dir}/${Head}_PB6_all_R1.pat.${body} |grep -v chrX | sort-bed --max-mem 10G - |gzip >${Head}_pat_temp.gz  &
done
wait

# Perform BinomialTest one strand at a time, using pooled MB6 and PB6 reads 
for Head in BN HT  SK  SP  KD  LV  GI  ST
do 
    MAT_READ_BED=${Head}_mat_temp.gz
    PAT_READ_BED=${Head}_pat_temp.gz
    BinomialTest_TSS ${Head}_allReads_TSS.bed ${Head}_allReads_TSS ${MAT_READ_BED} ${PAT_READ_BED} &
done

wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do 
j=${Head}_allReads_TSS
wc -l ${j}_binomtest_Rfdr0.1.bed ${j}_binomtest_Rfdr1.bed
done
  #  2423 BN_allReads_TSS_binomtest_Rfdr0.1.bed
  # 76850 BN_allReads_TSS_binomtest_Rfdr1.bed
  # 79273 total
  #  1291 HT_allReads_TSS_binomtest_Rfdr0.1.bed
  # 41606 HT_allReads_TSS_binomtest_Rfdr1.bed
  # 42897 total
  #  1109 SK_allReads_TSS_binomtest_Rfdr0.1.bed
  # 31708 SK_allReads_TSS_binomtest_Rfdr1.bed
  # 32817 total
  #  3426 SP_allReads_TSS_binomtest_Rfdr0.1.bed
  # 52414 SP_allReads_TSS_binomtest_Rfdr1.bed
  # 55840 total
  #  2550 KD_allReads_TSS_binomtest_Rfdr0.1.bed
  # 38101 KD_allReads_TSS_binomtest_Rfdr1.bed
  # 40651 total
  #  5793 LV_allReads_TSS_binomtest_Rfdr0.1.bed
  # 89980 LV_allReads_TSS_binomtest_Rfdr1.bed
  # 95773 total
  #  1753 GI_allReads_TSS_binomtest_Rfdr0.1.bed
  # 40308 GI_allReads_TSS_binomtest_Rfdr1.bed
  # 42061 total
  #  1454 ST_allReads_TSS_binomtest_Rfdr0.1.bed
  # 33340 ST_allReads_TSS_binomtest_Rfdr1.bed
  # 34794 total

#for Head in BN HT  SK  SP  KD  LV  GI  ST
#do 
#j=${Head}_allReads_TSS_maxTSNs_SNPs20bp
#cat ${j}_binomtest.bed                 | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_IGV.bed
#cat ${j}_binomtest_interestingHets.bed | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_interestingHets_IGV.bed
#cat ${j}_binomtest_Rfdr0.1.bed         | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_Rfdr0.1_IGV.bed  
#cat ${j}_binomtest_Rfdr1.bed           | awk 'BEGIN{OFS="\t"} (NR>1 && $11+0 >0.9){print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_Rfdr0.9_IGV.bed 
#cat ${j}_binomtest_Rfdr0.2.bed         | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_Rfdr0.2_IGV.bed  
#done


cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/maxTSN_TSS_TID_combine_analysis_MultiBaseRunOn
ln -s ../TSS_BinomialTest_MultiBaseRunOn/*_allReads_TSS_binomtest_Rfdr1.bed .
ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_snps_CAST.subsample.bam/P.CAST.EiJ_M.C57BL.6J_paternal_all.fa .
ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_snps_CAST.subsample.bam/P.CAST.EiJ_M.C57BL.6J_maternal_all.fa .
ln -s ../TSS_KStest_MultiBaseRunOn/*_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_pValue_maxTSNs.bed .

# one TSS might contain more than 1 maxTSN
for Head in BN HT  SK  SP  KD  LV  GI  ST
do 
bedtools intersect -wo -s -a ${Head}_${studyBed}_5mat5pat_uniq_maskedVSunmasked_pValue_maxTSNs.bed \
-b <(cat ${Head}_allReads_TSS_binomtest_Rfdr1.bed| awk '{OFS="\t"} NR>1 {print $1, $2, $3, $4, $11, $10}')\
| awk '{OFS="\t"} ($2==$19 && $3==$20) {print $1, $2, $3, $4, $5, $6, $7,$8,$9,$10, $21,$22, $11,$12,$13}' |sort-bed - \
> ${Head}_${studyBed}_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs.bed &

# col 1-12  TSS 
# col 4 TSN Counts
# col 5 total (?) reads count
# col 6 strand
# col 7-10 maksed pvalue, masked fdr, unmasked p value, unmasked fdr
# col 11, P is Cast, M is B6. P,53,87 = Cast has more reads (might not be significant) , 53 B6 reads, 87 Cast reads
# col 12 Binomial test fdr
# col 13-15 maxTNS location, 

done

wait

GC_content_TSS(){
 Head=$1
 j=$2
 d=$3
 step=$4

## use +- d bp to identify seq
# get the bed geions  
# get the sequence from fasta
# -s  Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.

if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($13,4)p, $14-d, $15+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($13,4)p, $14-d, $15+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_mat_patSeq.bed 
 cat ${j}_+-${d}_mat_patSeq.bed  | awk '{OFS="\t"} (substr($11,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($11,1,1),$11, $12, $13, $14, $16, $17} 
 (substr($11,1,1)=="P") {print  $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($11,1,1), $11, $12, $13, $14, $17, $16}' > ${j}_+-${d}_High_LowAlleleSeq.bed 
# ignore S
# (substr($11,1,1)=="S") {print $1,$2,$3,$4,$5, $6, $7, $8, $9, $10, substr($11,1,1),$11, $12, $13, $14, $16, $17} 

# col 1-12  TSS 
# col 4 TSN Counts
# col 5 total (?) reads count
# col 6 strand
# col 7-10 maksed pvalue, masked fdr, unmasked p value, unmasked fdr
# col 11, P is Cast, M is B6. P,53,87 = Cast has more reads (might not be significant) , 53 B6 reads, 87 Cast reads
# col 12 Binomial test fdr
# col 13-15 maxTNS location, 
# col 16 High Allele Seq (Higher reads count in TSS)
# col 17 Low Allele Seq

fi

  R --vanilla --slave --args ${j}_+-${d}_High_LowAlleleSeq.bed  ${Head} ${step} ${d} < getGC_content_HighLowAllele.R 
}


for Head in BN LV 
do
d=35
step=5
GC_content_TSS ${Head} ${Head}_${studyBed}_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs $d ${step} &
done

# panel37_BN_deltaGC_exclude_CA_step=5_d=35.pdf , panel38_BN_deltaGC_exclude_CA_step=5_d=35-SI.pdf
Rscript GC_content_SeqLogo_manuscript_figure.R

