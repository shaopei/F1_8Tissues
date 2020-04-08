
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
   cat ${Head}_allReads_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8 >=b) {print $1, $2+$7-1, $2+$7 , $7c$6C$8, $5, $6}' |sort-bed - > ${Head}_allReads_TSN${b}+_IGV.bed &
   cat ${Head}_allReads_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8 >=b) {print $1, $2+$7-1, $2+$7 , $7, $8, $6}' |sort-bed - > ${Head}_allReads_TSN_pos_readcount${b}+_strand.bed &
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
wait

## identify maxTSNs with EACH TSS
# use mapped reads from AlleleDB (bowtie), including both mat, pat and identical reads

# identify the abundance of PolII at each position of TSS
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
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

## use +1 10bp to identify seqlogo
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
# Identify maxTSNs with SNPs within 20bp, downstream only
unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered
cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }' |sort-bed - |gzip > unfiltered_snp.sorted.bed.gz

# examine the distance distribution of closest donwstream SNP to maxTSNs
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
 bedtools closest -iu -D a -a <(sort-bed ${Head}_allReads_TSS_maxTSNs.bed) -b <(zcat unfiltered_snp.sorted.bed.gz) > ${Head}_allReads_TSS_maxTSNs_DistToDownstreamSNP.bed &
done


# examine the distance distribution of closest donwstream SNP to maxTSNs
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
 R --vanilla --slave --args $(pwd) ${Head}_allReads_TSS_maxTSNs_DistToDownstreamSNP.bed 11 100 < getHistFromCol.R &
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




# identify the abundance of PolII at each position of TSS
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
for Head in BN HT  SK  SP  KD  LV  GI  ST
do
  bedtools coverage -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz ) > ${Head}_allReads_TSS_PosAllReadCounts.bed &
done


