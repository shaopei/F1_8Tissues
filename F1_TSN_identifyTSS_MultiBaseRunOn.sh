
# make 1bp bed file using the 5 prime of the nascent RNA
wd=/workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn
cd ${wd}
for folder in allelicbias-PersonalGenome_P.CAST_M.B6-*_R1
  do #echo ${folder}
  PREFIX=`echo ${folder}|rev |cut -d - -f 1 |rev`
  #echo ${PREFIX}
  echo "bash make_map2ref_1bpbed_map5_MultiBaseRunOn_B6_Cast_F1.bsh ${PREFIX} ${wd}/allelicbias-PersonalGenome_P.CAST_M.B6-${PREFIX} > make_map2ref_1bpbed_map5_MultiBaseRunOn_${PREFIX}.log 2>&1 &"
done
mkdir map2ref_1bpbed_map5
ln -s ../allelicbias-PersonalGenome_P.CAST_M.B6-*_all_R1/*.map2ref.map5.1bp.sorted.bed.gz .
# HERE

# identify TSB within dREG sites
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn
studyBed=dREG
ln -s /workdir/sc2457/F1_Tissues/dREG/Browser/ .
ln -s /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .


## identify the abundance of PolII at each position of TID (dREG)
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0.bed &
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*_identical.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0_ide.bed &
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0_mat.bed &
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0_pat.bed &
done

## identify TSN
# only keep base with at least b all reads ($8 >=b)
b=5
wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
   cat ${Head}_allReads_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8-0 >=b) {print $1, $2+$7-1, $2+$7 , $7c$6C$8, $5, $6}' |sort-bed - > ${Head}_allReads_TSN${b}+_IGV.bed &
   cat ${Head}_allReads_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8-0 >=b) {print $1, $2+$7-1, $2+$7 , $7, $8, $6}' |sort-bed - > ${Head}_allReads_TSN_pos_readcount${b}+_strand.bed &
   # $4 is TSN position inside dREG, $5 read counts of the TSN, $6 strand of the TSN
done

# only keep base with at least 7 mat (B6) or pat(cast) reads ($8 >=b)
#b=2
#wait
#for Head in BN HT  SK  SP  KD  LV  GI  ST
#do
#   cat ${Head}_allReads_temp0_pat.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8 >=b) {print $1, $2+$7-1, $2+$7 , $7c$6C$8, $5, $6}' |sort-bed - > ${Head}_allReads_TSN${b}+_mat_IGV.bed &
#   cat ${Head}_allReads_temp0_pat.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8 >=b) {print $1, $2+$7-1, $2+$7 , $7, $8, $6}' |sort-bed - > ${Head}_allReads_TSN_pos_readcount${b}+_strand.bed &
   # $4 is TSN position inside dREG, $5 read counts of the TSN, $6 strand of the TSN
#done



## identify TSS
# Force strandedness -s
# merge TSN with gap <= 60bp -d 60
wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
   bedtools merge -s -d 60 -c 4,5,6 -o count,sum,distinct -i ${Head}_allReads_TSN_pos_readcount${b}+_strand.bed |grep -v chrY > ${Head}_allReads_TSS.bed &
   # $4 is number of TSN in the TSS, $5 sum of the read counts of the TSN (with at least 2 reads), $6 strand of the TSS
done
wait

## identify maxTSNs with EACH TSS
# use Proseq2.0, BWA mapping to mm10
# output: "chr" "chrStart"  "chrEnd"  "TSNCount(ofTSS)"  "ReadsCount(sumOfQualifiedTSNReadsCount)"   "Strand"   "map5.peaks.posotion"   "maxReadCountOftheTSN"
#for Head in BN HT  SK  SP  KD  LV  GI  ST
#do
#	rm ${Head}_allReads_TSS_maxTSNsCol7.bed
#    R --vanilla --slave --args $(pwd)  ${Head} _allReads_TSS bigWig/ _all_ < getMaxTSN_cbsudanko.R &
#done


## identify maxTSNs with EACH TSS
# use mapped reads from AlleleDB (bowtie), including both mat, pat and identical reads

# identify the abundance of PolII at each position of TSS
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1.bed &
  #bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*_identical.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1_ide.bed &
  #bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1_mat.bed &
  #bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1_pat.bed &
done
wait

for Head in BN HT  SK  SP  KD  LV  GI  ST
do
   python getMaxTSNs_frombedtools_coverage_strandSpecific.py ${Head}_allReads_TSStemp1.bed ${Head}_allReads_TSS_maxTSNsCol7_minusStrandSameDirection.bed &
done
wait


for Head in BN HT  SK  SP  KD  LV  GI  ST
do
   #cat ${Head}_allReads_TSS_maxTSNsCol7.bed | awk '{OFS="\t"} ($6=="+"){print $1, $2+$7-1, $2+$7, $8, "111", $6} ($6=="-"){print $1, $3-$7, $3-$7+1, $8, "111", $6}' >  ${Head}_allReads_TSS_maxTSNs.bed &
   cat ${Head}_allReads_TSS_maxTSNsCol7_minusStrandSameDirection.bed | awk '{OFS="\t"} {print $1, $2+$7-1, $2+$7, $8, "111", $6}' >  ${Head}_allReads_TSS_maxTSNs.bed &

done

## use +- 10bp to identify seqlogo
wait
d=10
# get the bed geions 
# col4 is ReadCountOftheTSN
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
   cat ${Head}_allReads_TSS_maxTSNs.bed | awk -v d=$d '{OFS="\t"} {print $1, $2-d, $3+d, $4, $5, $6}' >  ${Head}_allReads_TSS_maxTSNs+-${d}.bed &
done

wait
# get the sequence from fasta
# -s	Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
 bedtools getfasta -s -fi mm10.fa -bed ${Head}_allReads_TSS_maxTSNs+-${d}.bed  | grep -v \> > ${Head}_allReads_TSS_maxTSNs+-${d}.txt &
done
wait
# the the seqlogo
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
   R --vanilla --slave --args $(pwd) ${Head}_allReads_TSS_maxTSNs+-10.txt ${Head}_allReads_TSS_maxTSNs+-10_SeqLogo.pdf < getSeqLogo.R &
done



### identify maxTSN tagged with SNPS
unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered
cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }' |sort-bed - |gzip > unfiltered_snp.sorted.bed.gz

# examine the distance distribution of closest donwstream SNP to maxTSNs
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
 bedtools closest -iu -D a -a <(sort-bed ${Head}_allReads_TSS_maxTSNs.bed) -b <(zcat unfiltered_snp.sorted.bed.gz) > ${Head}_allReads_TSS_maxTSNs_DistToDownstreamSNP.bed &
done
wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
 R --vanilla --slave --args $(pwd) ${Head}_allReads_TSS_maxTSNs_DistToDownstreamSNP.bed 11 100 < getHistFromCol.R &
done

# Identify maxTSNs with SNPs within 20bp, downstream only
l=20
for Head in BN HT  SK  SP  KD  LV  GI  ST
 do 
cat ${Head}_allReads_TSS_maxTSNs_DistToDownstreamSNP.bed | awk -v l=$l '{OFS="\t"} ($11-0==$11 && $11-0 <= l && $11 != -1){print $0}' | cut -f 1-6 > ${Head}_allReads_TSS_maxTSNs_SNPs${l}bp.bed &
done



### what is the distribution of allelic read length?
cd /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/mappedReadLength
for f in *_AMBremoved_sorted_specific.bed.gz
do
zcat $f | awk '{OFS="\t"} {print $0, $3-$2}' > ${f}_temp2  &
done
wait
for f in *_AMBremoved_sorted_specific.bed.gz
do
R --vanilla --slave --args $(pwd) ${f}_temp2 7 100 < getHistFromCol.R &
done

for f in *_all_R1.mat.bowtie.gz_AMBremoved_sorted_identical.bed.gz
do
zcat $f | awk '{OFS="\t"} {print $0, $3-$2}' > ${f}_temp2  &
done
wait
for f in *_all_R1.mat.bowtie.gz_AMBremoved_sorted_identical.bed.gz
do
R --vanilla --slave --args $(pwd) ${f}_temp2 7 100 < getHistFromCol.R &
done


# use allelic reads that are at least lbp (>=l)
l=20
wd=/workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn
cd ${wd}
echo "keep read length >= ${l}"
for folder in allelicbias-PersonalGenome_P.CAST_M.B6-*_R1
  do #echo ${folder}
  PREFIX=`echo ${folder}|rev |cut -d - -f 1 |rev`
  #echo ${PREFIX}
  echo "bash make_map2ref_1bpbed_map5_ReadLengthFiltered.bsh ${PREFIX} ${wd}/allelicbias-PersonalGenome_P.CAST_M.B6-${PREFIX} ${l} > make_map2ref_1bpbed_map5_ReadLengthFiltered_${PREFIX}.log 2>&1 &"
done

mkdir map2ref_1bpbed_map5_FromMappedReadAtLeast20bpLong
cd ${wd}/map2ref_1bpbed_map5_FromMappedReadAtLeast20bpLong
ln -s ../allelicbias-PersonalGenome_P.CAST_M.B6-*/*map2ref.map5.1bp.R20bp.sorted.bed.gz .


cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn
ln -s ${wd}/map2ref_1bpbed_map5_FromMappedReadAtLeast20bpLong .

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
  python ${PL}/BinomialTestFor_merged_cov.bed.py ${j}.merged_cov.bed ${j}_binomtest.bed
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 0.1 ${j}_binomtest_Rfdr0.1.bed < getCorrectedPValue.R &
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 0.2 ${j}_binomtest_Rfdr0.2.bed < getCorrectedPValue.R &
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed 9 1 ${j}_binomtest_Rfdr1.bed < getCorrectedPValue.R &
  python ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt 
  awk 'NR==1 { print $0 } NR>1 && ($9+0) <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_binomtest_interestingHets.bed
  mv ${j}.mat_cov.bed ${j}.pat_cov.bed ${j}.merged_cov.bed toremove
}

# Pool MB6 and PB6 reads together. mat is B6 reads, pat is Cast reads
# remove chrX
bed_dir=map2ref_1bpbed_map5_FromMappedReadAtLeast20bpLong
body=bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.R20bp.sorted.bed.gz
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
    BinomialTest_TSN ${Head}_allReads_TSS_maxTSNs_SNPs20bp.bed ${Head}_allReads_TSS_maxTSNs_SNPs20bp ${MAT_READ_BED} ${PAT_READ_BED} &
done

wait
for Head in BN HT  SK  SP  KD  LV  GI  ST
do 
j=${Head}_allReads_TSS_maxTSNs_SNPs20bp
wc -l ${j}_binomtest_Rfdr0.1.bed ${j}_binomtest_Rfdr0.2.bed ${j}_binomtest_interestingHets.bed
done

for Head in BN HT  SK  SP  KD  LV  GI  ST
do 
j=${Head}_allReads_TSS_maxTSNs_SNPs20bp
cat ${j}_binomtest.bed                 | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_IGV.bed
cat ${j}_binomtest_interestingHets.bed | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_interestingHets_IGV.bed
cat ${j}_binomtest_Rfdr0.1.bed         | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_Rfdr0.1_IGV.bed  
cat ${j}_binomtest_Rfdr1.bed           | awk 'BEGIN{OFS="\t"} (NR>1 && $11+0 >0.9){print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_Rfdr0.9_IGV.bed 
cat ${j}_binomtest_Rfdr0.2.bed         | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, "111", $10}'  > ${j}_binomtest_Rfdr0.2_IGV.bed  
done

# examine the distribution of SNPs near maxTSN.
# Compare as.maxTSN (R, fdr<=0.1) vs non.as.maxTSN (R, fdr>= 0.9)
Rscript getSNPsAbundance.R

### seperate allelic different maxTSN into two groups: one inside allelic shape different TSS, one ouside
# use ${j}_binomtest_Rfdr0.1.bed  
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/maxTSN_TSS_TID_combine_analysis_MultiBaseRunOn
ln -s ../TSS_KStest_MultiBaseRunOn/*_allReads_TSS_5mat5pat_uniq_fdr0.1.bed .
ln -s ../identifyTSS_maxTSNs_MultiBaseRunOn/*_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_IGV.bed .
ln -s ../TID.dREG_KStest_MultiBaseRunOn/*_dREG_5mat5pat_uniq_fdr0.1.bed .

# examine if as.maxTSN is inside asTSS or asTID, strand specific
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -s -a ${Head}_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_IGV.bed -b ${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed > ${Head}_maxTSN_AS.TSS_temp &
  bedtools coverage -sorted -s -a ${Head}_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_IGV.bed -b ${Head}_dREG_5mat5pat_uniq_fdr0.1.bed > ${Head}_maxTSN_AS.TID_temp &
done

for Head in BN HT  SK  SP  KD  LV  GI  ST
do
paste ${Head}_maxTSN_AS.TSS_temp  ${Head}_maxTSN_AS.TID_temp | awk '{OFS="\t"} ($2==$12) {print $1,$2,$3,$4,$5,$6, $7, $17}' > ${Head}_maxTSN_AsTSS_AsTID.bed
mv ${Head}_maxTSN_AS.TSS_temp  ${Head}_maxTSN_AS.TID_temp toremove/.
done

echo -e "Organ\tmaxTSN\tmaxTSN_asTSS_only\tmaxTSN_asTID_only\tmaxTSN_asTSS_AND_asTID" > maxTSN_AsTSS_AsTID_summary.txt
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  # maxTSN, maxTSN_asTSS_only, maxTSN_asTID_only, maxTSN_asTSS_AND_asTID
  echo ${Head} > temp
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | wc -l  >> temp
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | awk '{OFS="\t"} ($7==1 && $8==0){print $0}' |wc -l >> temp
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | awk '{OFS="\t"} ($8==1 && $7==0){print $0}' |wc -l >> temp
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | awk '{OFS="\t"} ($7+$8==2){print $0}' |wc -l >> temp
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | awk '{OFS="\t"} ($7+$8==1){print $0}' |wc -l >> temp

  paste temp -s >> maxTSN_AsTSS_AsTID_summary.txt
done

# uniq TID (not just allelic different one) that tagged with as.maxTSN
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -sorted -s -a ${Head}_dREG_5mat5pat_uniq.bed -b ${Head}_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_IGV.bed \
  | awk '{OFS="\t"} ($7>=1){print $0}' | cut -f 1-6 |sort-bed - |uniq
done
# the subset of uniq TID with 2+ TSS(that tagged with SNPs)
# the asTID that


# check some exmaples at IGV
Head=BN
#maxTSN_with_asTSS_only (no asTID)
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | awk '{OFS="\t";m=":"; d="-"} ($7==1 && $8==0){print $0, $1m$2d$3}'
  #maxTSN_with_asTID_only (no asTSS)
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | awk '{OFS="\t";m=":"; d="-"}  ($8==1 && $7==0){print $0, $1m$2d$3}'
  #maxTSN with Both asTID and asTSS
  cat ${Head}_maxTSN_AsTSS_AsTID.bed | awk '{OFS="\t";m=":"; d="-"}  ($8==1 && $7==1){print $0, $1m$2d$3}'

# asTID with at least 2 asTSS (not check maxTSN, the only 1 of the output contains maxTSN in brain)
echo -e "Organ\tasTIDwith2+asTSS_Count" > AsTID_w2+asTSS_summary.txt
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  echo ${Head} > temp
bedtools coverage -sorted -a ${Head}_dREG_5mat5pat_uniq_fdr0.1.bed   -b ${Head}_allReads_TSS_5mat5pat_uniq_fdr0.1.bed | awk '{OFS="\t";m=":"; d="-"} ($7 >= 2) {print $0, $1m$2d$3}' |cut -f 1-3 |uniq |wc -l >> temp
paste temp -s >> AsTID_w2+asTSS_summary.txt
done


### calculate Tm around maxTSN. High allele and low allele
# generate maternal and paternal genome to get the sequence 
cd /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_snps_CAST.subsample.bam
rm P.CAST.EiJ_M.C57BL.6J_paternal_all.fa
for f in *_P.CAST.EiJ_M.C57BL.6J_paternal.fa
do cat $f >> P.CAST.EiJ_M.C57BL.6J_paternal_all.fa
done

cd /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_snps_CAST.subsample.bam
rm P.CAST.EiJ_M.C57BL.6J_maternal_all.fa
for f in *_P.CAST.EiJ_M.C57BL.6J_maternal.fa
do cat $f >> P.CAST.EiJ_M.C57BL.6J_maternal_all.fa
done

cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_MultiBaseRunOn
ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_snps_CAST.subsample.bam/P.CAST.EiJ_M.C57BL.6J_paternal_all.fa .
ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_snps_CAST.subsample.bam/P.CAST.EiJ_M.C57BL.6J_maternal_all.fa .
 
## use +- d bp to identify seq
wait
d=20
# get the bed geions 
# col4 and col5 is win_pareant,mat_read_count,pat_read_count (col5 is S if not significantly allelic biased)
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  j=${Head}_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1
   cat ${j}.bed | awk -v d=$d '{OFS="\t"} NR>1 {print $1, $2-d, $3+d, $4, $5, $10}' >  ${j}_+-${d}.bed &
done

wait
# get the sequence from fasta
# -s  Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.


for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  j=${Head}_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1
 #bedtools getfasta -s -fi mm10.fa -bed ${j}_+-${d}.bed  | grep -v \> > ${j}_mm10.txt &
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}_+-${d}.bed |awk '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2,$3,$4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}_+-${d}.bed |awk '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2,$3,$4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
done
wait

# examine if the result is the same as using mm10
# diff -i  BN_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_mm10.txt BN_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_P.CAST.EiJ_M.C57BL.6J_maternal.txt
# YES!

for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  j=${Head}_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1
  paste ${j}_+-${d}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_mat_patSeq.bed 
  cat ${j}_+-${d}_mat_patSeq.bed  | awk '{OFS="\t"} (substr($4,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8} (substr($4,1,1)=="P") {print $1,$2,$3,$4,$5, $6, $8, $7}' > ${j}_+-${d}_High_LowAlleleSeq.bed 
done

Rscript getGC_content_HighLowAllele.R




