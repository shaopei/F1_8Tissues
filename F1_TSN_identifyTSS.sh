# identify TSB within dREG sites
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS

studyBed=dREG
ln -s /workdir/sc2457/F1_Tissues/dREG/Browser/ .
ln -s /workdir/sc2457/F1_Tissues/SingleBaseRunOn/map2ref_1bpbed_map5 .
ln -s /workdir/sc2457/F1_Tissues/SingleBaseRunOn/map2ref_1bpbed_map5 map2ref_1bpbed


## identify the abundance of PolII at each position
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
for Head in HT KD SK
do
  bedtools coverage -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz ) > ${Head}_allReads_temp0.bed &
  #bedtools coverage -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*_identical.map2ref.map5.1bp.sorted.bed.gz ) > ${Head}_allReads_temp0_ide.bed &
  #bedtools coverage -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz ) > ${Head}_allReads_temp0_mat.bed &
  #bedtools coverage -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz ) > ${Head}_allReads_temp0_pat.bed &
done

# only keep base with at least 2 all reads ($8 >1)
wait
for Head in HT KD SK
do
   cat ${Head}_allReads_temp0.bed | awk 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8 >1) {print $1, $2+$7-1, $2+$7 , $7c$6C$8, $5, $6}' |sort-bed - > ${Head}_allReads_TSN_IGV.bed &
   cat ${Head}_allReads_temp0.bed | awk 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8 >1) {print $1, $2+$7-1, $2+$7 , $7, $8, $6}' |sort-bed - > ${Head}_allReads_TSN_pos_readcount_strand.bed &
   # $4 is TSN position inside dREG, $5 read counts of the TSN, $6 strand of the TSN
done

## identify TSS
# Force strandedness -s
# merge TSN with gap <= 60bp -d 60
wait
for Head in HT KD SK
do
   bedtools merge -s -d 60 -c 4,5,6 -o count,sum,distinct -i ${Head}_allReads_TSN_pos_readcount_strand.bed > ${Head}_allReads_TSS.bed &
   # $4 is number of TSN in the TSS, $5 sum of the read counts of the TSN (with at least 2 reads), $6 strand of the TSS
done

## identify maxTSN with EACH TSS
# use Proseq2.0, BWA mapping to mm10
# output: "chr" "chrStart"  "chrEnd"  "TSNCount(ofTSS)"  "ReadsCount(sumOfQualifiedTSNReadsCount)"   "Strand"   "map5.peaks.posotion"   "maxReadCountOftheTSN"
for Head in HT KD SK
do
   R --vanilla --slave --args $(pwd) map5_bw/ ${Head} _allReads_TSS < getMaxTSN_cbsudanko.R &
done

## use +1 10bp to identify seqlogo
wait
d=10
# get the bed geions 
# col4 is ReadCountOftheTSN
for Head in HT KD SK
do
   cat ${Head}_allReads_TSS_maxTSNs.bed | awk -v d=$d '{OFS="\t"} {print $1, $2+$7-1-d, $2+$7+d, $8, "111", $6}' >  ${Head}_allReads_TSS_maxTSNs+-${d}.bed &
done
wait
# get the sequence from fasta
# -s	Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.
for Head in HT KD SK
do
 bedtools getfasta -s -fi mm10.fa -bed ${Head}_allReads_TSS_maxTSNs+-${d}.bed  | grep -v \> > ${Head}_allReads_TSS_maxTSNs+-${d}.txt &
done
# the the seqlogo
for Head in HT KD SK
do
   R --vanilla --slave --args $(pwd) ${Head}_allReads_TSS_maxTSNs+-10.txt ${Head}_allReads_TSS_maxTSNs+-10_SeqLogo.pdf < getSeqLogo.R &
done












# identify the abundance of PolII at each position of TSS
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
for Head in HT KD SK
do
  bedtools coverage -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz ) > ${Head}_allReads_TSS_PosAllReadCounts.bed &
done

