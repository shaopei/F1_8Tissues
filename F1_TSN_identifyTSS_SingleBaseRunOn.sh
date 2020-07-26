### used for single base runon ChRO-seq
# identify TSB within dREG sites
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_SingleBaseRunOn

ln -s /workdir/sc2457/F1_Tissues/dREG/Browser/ .
ln -s /workdir/sc2457/F1_Tissues/SingleBaseRunOn/map2ref_1bpbed_map5 .
ln -s /workdir/sc2457/F1_Tissues/SingleBaseRunOn/map2ref_1bpbed_map5 map2ref_1bpbed

## identify the abundance of PolII at each position of TID (dREG)
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
for Head in HT KD SK
do
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0.bed &
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*_identical.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0_ide.bed &
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0_mat.bed &
  bedtools coverage -sorted -d -s -a <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz| awk 'BEGIN{OFS="\t"} {print $0, ".", "+"} {print $0, ".", "-"}') -b <(zcat map2ref_1bpbed_map5/${Head}*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_temp0_pat.bed &
done
# HERE
## identify TSN
# only keep base with at least b all reads ($8 >=b)
b=2
wait
for Head in HT KD SK
do
   cat ${Head}_allReads_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8-0 >=b) {print $1, $2+$7-1, $2+$7 , $7c$6C$8, $5, $6}' |sort-bed - > ${Head}_allReads_TSN${b}+_IGV.bed &
   cat ${Head}_allReads_temp0.bed | awk -v b=$b 'BEGIN{OFS="\t"; c=","; C=",reads:"} ($8-0 >=b) {print $1, $2+$7-1, $2+$7 , $7, $8, $6}' |sort-bed - > ${Head}_allReads_TSN_pos_readcount${b}+_strand.bed &
   # $4 is TSN position inside dREG, $5 read counts of the TSN, $6 strand of the TSN
done

## identify TSS
# Force strandedness -s
# merge TSN with gap <= 60bp -d 60
wait
for Head in HT KD SK
do
   bedtools merge -s -d 60 -c 4,5,6 -o count,sum,distinct -i ${Head}_allReads_TSN_pos_readcount${b}+_strand.bed |grep -v chrY > ${Head}_allReads_TSS.bed &
   # $4 is number of TSN in the TSS, $5 sum of the read counts of the TSN (with at least b reads), $6 strand of the TSS
done
wait


## identify maxTSNs with EACH TSS
# use mapped reads from AlleleDB (bowtie), including both mat, pat and identical reads

# identify the abundance of PolII at each position of TSS
# strand specific
# use all reads (not just allelic reads), 
# use mapping from bowtie mapping, pat is liftovered 
wait
for Head in HT KD SK
do
  bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1.bed &
  bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*_identical.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1_ide.bed &
  bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1_mat.bed &
  bedtools coverage -sorted -d -s -a ${Head}_allReads_TSS.bed -b <(zcat map2ref_1bpbed_map5/${Head}*.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz |sort-bed --max-mem 10G -) > ${Head}_allReads_TSStemp1_pat.bed &
done
wait


# paste the all reads  and alleleic reads count into a table
# only keep base with at least b all reads ($8 >= k)
k=5
for Head in HT KD SK
do
  for allele in mat pat
  do
    paste ${Head}_allReads_TSStemp1.bed ${Head}_allReads_TSStemp1_${allele}.bed \
    | awk -v k=$k 'BEGIN{OFS="\t"} ($2==$10 && $3==$11 && $8+0 >= k ) {print $0}' |cut -f 9- > ${Head}_${allele}_temp2.bed &
  done
done


# report more than one maxTSN if multiple TSN share the same max read count
for Head in HT KD SK
do
   python2 getMaxTSNs_frombedtools_coverage_strandSpecific.py ${Head}_allReads_TSStemp1.bed ${Head}_allReads_TSS_maxTSNsCol7_minusStrandSameDirection.bed &
done

wait

for Head in HT KD SK
do
   #cat ${Head}_allReads_TSS_maxTSNsCol7.bed | awk '{OFS="\t"} ($6=="+"){print $1, $2+$7-1, $2+$7, $8, "111", $6} ($6=="-"){print $1, $3-$7, $3-$7+1, $8, "111", $6}' >  ${Head}_allReads_TSS_maxTSNs.bed &
   cat ${Head}_allReads_TSS_maxTSNsCol7_minusStrandSameDirection.bed | awk '{OFS="\t"} {print $1, $2+$7-1, $2+$7, $8, "111", $6}' >  ${Head}_allReads_TSS_maxTSNs.bed &
done

## use +- 10bp to identify seqlogo
wait
d=10
# get the bed geions 
# col4 is ReadCountOftheTSN
for Head in HT KD SK
do
   cat ${Head}_allReads_TSS_maxTSNs.bed | awk -v d=$d '{OFS="\t"} {print $1, $2-d, $3+d, $4, $5, $6}' >  ${Head}_allReads_TSS_maxTSNs+-${d}.bed &
done

wait
# get the sequence from fasta
# -s	Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.
for Head in HT KD SK
do
 bedtools getfasta -s -fi mm10.fa -bed ${Head}_allReads_TSS_maxTSNs+-${d}.bed  | grep -v \> > ${Head}_allReads_TSS_maxTSNs+-${d}.txt &
done
wait
# the the seqlogo
for Head in HT KD SK
do
   R --vanilla --slave --args $(pwd) ${Head}_allReads_TSS_maxTSNs+-10.txt ${Head}_allReads_TSS_maxTSNs+-10_SeqLogo.pdf < getSeqLogo.R &
done



# maxTSN of allelic reads
# report more than one maxTSN if multiple TSN share the same max read count
for Head in HT KD SK
do
  for allele in mat pat
  do
   python2 getMaxTSNs_frombedtools_coverage_strandSpecific.py ${Head}_${allele}_temp2.bed ${Head}_${allele}Reads_TSS_maxTSNsCol7_minusStrandSameDirection.bed &
 done
done
wait


# only keep allelic maxTSN with at least 5 allelic reads
for Head in HT KD SK
do
  for allele in mat pat
  do
    cat ${Head}_${allele}Reads_TSS_maxTSNsCol7_minusStrandSameDirection.bed  | awk '{OFS="\t"} ($8+0>5){print $1, $2+$7-1, $2+$7, $8, "111", $6}' \
     >  ${Head}_${allele}Reads_TSS_maxTSNs.bed  # col4 is read count
        mv ${Head}_${allele}Reads_TSS_maxTSNsCol7_minusStrandSameDirection.bed toremove/.
  done
done

# identify allelic maxTSN that share the same site (overlap)
for Head in HT KD SK
do intersectBed -wo -s -a ${Head}_matReads_TSS_maxTSNs.bed -b ${Head}_patReads_TSS_maxTSNs.bed > ${Head}_matReads_patReads_TSS_maxTSNs.bed &
done

[sc2457@cbsudanko identifyTSS_SingleBaseRunOn]$ wc -l *atReads_TSS_maxTSNs.bed
   5202 HT_matReads_TSS_maxTSNs.bed
   1976 HT_matReads_patReads_TSS_maxTSNs.bed
   5055 HT_patReads_TSS_maxTSNs.bed
   4474 KD_matReads_TSS_maxTSNs.bed
   1865 KD_matReads_patReads_TSS_maxTSNs.bed
   4464 KD_patReads_TSS_maxTSNs.bed
   4995 SK_matReads_TSS_maxTSNs.bed
   1933 SK_matReads_patReads_TSS_maxTSNs.bed
   4909 SK_patReads_TSS_maxTSNs.bed
  34873 tota

for Head in HT KD SK
do cat ${Head}_matReads_patReads_TSS_maxTSNs.bed | awk '{OFS="\t"; c=","} ($4/$10<2 && $4/$10>0.5){print $1,$2,$3,$4c$10, $5, $6}' > ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
done

[sc2457@cbsudanko identifyTSS_SingleBaseRunOn]$ wc -l *_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  1448 HT_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  1395 KD_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  1461 SK_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  4304 total

# get the read
for Head in HT KD SK
do 
  for allele in mat pat
  do
    intersectBed -wb -a ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed -b <(zcat map2ref_1bpbed_map5/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz) \
    >  ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_reads.bed&
  done
done





