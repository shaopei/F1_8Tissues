### used for single base runon ChRO-seq

# insert(read) length distribution
for Head in HT KD SK
do
  for allele in mat pat
  do
    zcat  map2ref_1bpbed_map5/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz \
    | awk 'BEGIN{OFS="\t"}  {split($4,a,","); print length(a[2])}' > ${Head}_${allele}_specific_read_length.txt &
  done
  zcat  map2ref_1bpbed_map5/${Head}_PB6_*_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.map5.1bp.sorted.bed.gz \
  | awk 'BEGIN{OFS="\t"}  {split($4,a,","); print length(a[2])}' > ${Head}_mat_identical_read_length.txt &
done


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
# sup figure 4A
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
# paste the all reads  and alleleic reads count into a table
# only keep base with at least k all reads ($8 >= k)
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
   1976 HT_matReads_patReads_TSS_maxTSNs.bed  # shared maxTSNs between B6 and CAST allele
   5055 HT_patReads_TSS_maxTSNs.bed
   4474 KD_matReads_TSS_maxTSNs.bed
   1865 KD_matReads_patReads_TSS_maxTSNs.bed
   4464 KD_patReads_TSS_maxTSNs.bed
   4995 SK_matReads_TSS_maxTSNs.bed
   1933 SK_matReads_patReads_TSS_maxTSNs.bed
   4909 SK_patReads_TSS_maxTSNs.bed
  34873 total

for Head in HT KD SK
do cat ${Head}_matReads_patReads_TSS_maxTSNs.bed | awk '{OFS="\t"; c=","} ($4/$10<2 && $4/$10>0.5){print $1,$2,$3,$4c$10, $5, $6}' > ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
done

[sc2457@cbsudanko identifyTSS_SingleBaseRunOn]$ wc -l *_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  1448 HT_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  1395 KD_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  1461 SK_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed
  4304 total

# get the read that with 5 prime end mapped at the maxTSN of interests
for Head in HT KD SK
do 
  for allele in mat pat
  do
    intersectBed -s -wb -a ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed -b <(zcat map2ref_1bpbed_map5/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz) \
    >  ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_${allele}reads.bed &
  done
done

#head ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_${allele}reads.bed
#chr1    6214407 6214408 11,8    111     -       chr1    6214407 6214408 10068178,GCTCACTCGTCAGCCTCCGGTTCCCCT    111     -
#chr1    6214407 6214408 11,8    111     -       chr1    6214407 6214408 11224680,ACTCGTCAGCCTCCGGTTCCCCT        111     -
#chr1    6214407 6214408 11,8    111     -       chr1    6214407 6214408 12187797,CTCGTCAGCCTCCGGTTCCCCT 111     -
#chr1    6214407 6214408 11,8    111     -       chr1    6214407 6214408 12704199,CACTCGTCAGCCTCCGGTTCCCCT       111     -
#chr1    6214407 6214408 11,8    111     -       chr1    6214407 6214408 2446188,GGTGCTCACTCGTCAGCCTCCGGTTCCCCT  111     -

for Head in HT KD SK
do 
  for allele in mat pat
  do 
    python2 Generate_vector_input_for_KStest_fromReadLength.py \
    ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_${allele}reads.bed \
    ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_${allele}.ReadLength.bed 
  done
done

for Head in HT KD SK
do 
  R --vanilla --slave --args $(pwd) ${Head} matReads_patReads_TSS_maxTSNs_ratio0.5-2 < KStest_flexible_length_ReadLength.R &
done


# mapped location of pause
ln -s /workdir/sc2457/F1_Tissues/SingleBaseRunOn/map2ref_1bpbed_map3 .
cd map2ref_1bpbed_map3
for f in *AMBremoved_sorted_specific.map2ref.1bp.sorted.bed.gz
do j=`echo $f| rev|cut -d . -f 2-|rev` 
echo $j
  zcat $f > $j &
done
cd .. 

wait_a_second() {
  joblist=($(jobs -p))
    while (( ${#joblist[*]} >= 600 ))
      do
      sleep 1
      joblist=($(jobs -p))
  done
}

for Head in HT KD SK
do 
  for allele in mat pat
  do 
    rm ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_temp -r
    mkdir ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_temp
    while read  chr start end strand r
    do
  #echo $chr $start $end $strand $r
  #example of $r 10068178,GCTCACTCGTCAGCCTCCGGTTCCCCT
  #use $r to identify the reads from map2ref_1bpbed_map3 to indentify the mapped position of the reads at 3 prime end (pause)
  nice grep $r map2ref_1bpbed_map3/${Head}_PB6_F*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted.bed | cut -d ":" -f 2 \
  > ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_temp/${chr}_${start}_${end}_${strand}_$r &
  wait_a_second
    done < <(cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_${allele}reads.bed | cut -f 1-3,6,10)
  done
done

for Head in HT KD SK
do 
  for allele in mat pat
  do 
    rm ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads -r 
    mkdir ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads
    while read  chr start end strand
    do
  #echo $chr $start $end $strand
  # merge reads shared the same maxTSN into a file
  cat ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_temp/${chr}_${start}_${end}_${strand}_* >> ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads/${chr}_${start}_${end}_${strand}
done < <(cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_${allele}reads.bed | cut -f 1-3,6 | uniq)
done
done

for Head in HT KD SK
do 
  for allele in mat pat
  do 
    rm ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_map2ref_DistanceTomaxTSN.temp
    while read  chr start end  name1 name2 strand
    do
  # ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_map2refDistanceTomaxTSN/${chr}_${start}_${end}_${strand}
  # get the distance between map3(3 prime end of read) to maxTSN(-b, 5 prime end of reads)
  m=`bedtools closest -s -d -a <(sort-bed ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads/${chr}_${start}_${end}_${strand}) -b <(echo  -e "$chr\t$start\t$end\t$name1\t$name2\t$strand")  | cut -f 13 `
  if [[ "$m" == "" ]] ; then
    echo  -e "$chr\t$start\t$end\t$name1\t$name2\t$strand\tNA" >> ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_map2ref_DistanceTomaxTSN.temp
  else
  echo  -e "$chr\t$start\t$end\t$name1\t$name2\t$strand\t"$m >> ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_map2ref_DistanceTomaxTSN.temp
  fi
done < <(cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_${allele}reads.bed | cut -f 1-6 | uniq )

sed 's/ /,/g' ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_map2ref_DistanceTomaxTSN.temp > ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_map2ref_DistanceTomaxTSN.bed 
done
done

#KS test to identify maxTSNs that have allelic difference in pause shape
for Head in HT KD SK
do 
R --vanilla --slave --args $(pwd) ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_matreads_map2ref_DistanceTomaxTSN.bed  ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_patreads_map2ref_DistanceTomaxTSN.bed \
${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map3TomaxTSN_PValue.bed  < KStest_flexible_length_mat_pat_output.R &

done

# make map3 IGC track with only reads from the shared maxTSNs
#Head=KD
for allele in mat pat
  do 
cat ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads/* | sort-bed -  > ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads.bed &
done

#HERE
cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_SingleBaseRunOn

for Head in HT KD SK
do 
# make a file with both ReadLength and map3TomaxTSNDistance exclude any region contain NA (KS test input)
# NA is from the map3 position were lost during pat2ref liftover
intersectBed -wo -a <(paste ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_mat.ReadLength.bed ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_pat.ReadLength.bed  | cut -f 1-7,14) \
-b <(paste ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_matreads_map2ref_DistanceTomaxTSN.bed ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_patreads_map2ref_DistanceTomaxTSN.bed | cut -f 1-7,14) \
|cut -f 1-8,15,16 | grep -v NA > ${Head}_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed &
done

#wc -l *_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed
#   1446 HT_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed
#   1394 KD_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed
#   1461 SK_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed

# make a file with indel length with  ${Head}_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed
ln -s /workdir/sc2457/F1_Tissues/Pause_SingleBaseRunOn/pause_Indel/P.CAST_M.B6_F1hybrid.indels.bed .
for Head in HT KD SK
do 
bedtools closest -iu -D a -a <(sort -k1,1 -k2,2n ${Head}_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed) -b P.CAST_M.B6_F1hybrid.indels.bed > ${Head}_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed &
done

#wc -l *_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed
#   1446 HT_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed
#   1394 KD_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed
#   1461 SK_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed



# make a file with both ReadLength and map3TomaxTSNDistance (p value)
intersectBed -wo -a ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_ReadLengthPValue.bed -b ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map3TomaxTSN_PValue.bed \
|cut -f 1-8,15-17 > ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_ReadLengthPValue_map3TomaxTSNPValue.bed


cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_ReadLengthPValue_map3TomaxTSNPValue.bed | awk '{OFS="\t"; m=":"; d="-"} ($8+0>0.1 && $11+0<0.1){print $0, $1m$2d$3}' | sort -k 8

cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_ReadLengthPValue_map3TomaxTSNPValue.bed | awk '{OFS="\t"; m=":"; d="-"} ($8+0<0.1 && $11+0>0.1){print $0, $1m$2d$3}' | sort -k 11nr

cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_ReadLengthPValue_map3TomaxTSNPValue.bed | awk '{OFS="\t"; m=":"; d="-"} ($8+0<0.1 && $11+0<0.1){print $0, $1m$2d$3}' | sort -k 11nr



cd /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/identifyTSS_SingleBaseRunOn/Pasue_SNP_analysis
ln -s ../*_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed  .  # shared maxTSNs between B6 and CAST allele with allelic read ratio between 0.5-2

# identify maxPause using all reads from the maxTSNs
# get the read (including idnetical reads in mat)
for Head in HT KD SK
do 
  for allele in mat
  do
    intersectBed -s -wb -a ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2.bed -b <(zcat ../map2ref_1bpbed_map5/${Head}_PB6_*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_identical.map2ref.map5.1bp.sorted.bed.gz) \
    >  ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_IDEreads.bed &
  done
done

for Head in HT KD SK
do 
  for allele in mat pat
  do
ln -s ../${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_${allele}reads.bed .
  done
done


# get the map3 position for identical reads
allele=mat

wait_a_second() {
  joblist=($(jobs -p))
    while (( ${#joblist[*]} >= 600 ))
      do
      sleep 1
      joblist=($(jobs -p))
  done
}

for Head in HT KD SK
do 
    rm ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_temp -r
    mkdir ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_temp
    while read  chr start end strand r
    do
  #echo $chr $start $end $strand $r
  nice grep $r ../map2ref_1bpbed_map3/${Head}_PB6_F*_dedup_R1.${allele}.bowtie.gz_AMBremoved_sorted_identical.map2ref.1bp.sorted.bed | cut -d ":" -f 2 \
  > ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_temp/${chr}_${start}_${end}_${strand}_$r &
  wait_a_second
    done < <(cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_IDEreads.bed | cut -f 1-3,6,10)
  done

allele=mat
for Head in HT KD SK
do 
    rm ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads -r 
    mkdir ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads
    while read  chr start end strand
    do
  #echo $chr $start $end $strand
  cat ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_temp/${chr}_${start}_${end}_${strand}_* >> ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads/${chr}_${start}_${end}_${strand}
done < <(cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_IDEreads.bed | cut -f 1-3,6 | uniq)
done



for Head in HT KD SK
do 
    rm ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_map2ref_DistanceTomaxTSN.temp
    while read  chr start end  name1 name2 strand
    do
  m=`bedtools closest -s -d -a <(sort-bed ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads/${chr}_${start}_${end}_${strand}) -b <(echo  -e "$chr\t$start\t$end\t$name1\t$name2\t$strand")  | cut -f 13 `
  if [[ "$m" == "" ]] ; then
    echo  -e "$chr\t$start\t$end\t$name1\t$name2\t$strand\tNA" >> ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_map2ref_DistanceTomaxTSN.temp
  else
  echo  -e "$chr\t$start\t$end\t$name1\t$name2\t$strand\t"$m >> ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_map2ref_DistanceTomaxTSN.temp
  fi
done < <(cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_IDEreads.bed | cut -f 1-6 | uniq )

sed 's/ /,/g' ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_map2ref_DistanceTomaxTSN.temp > ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_map2ref_DistanceTomaxTSN.bed 
done


for Head in HT KD SK
do 
  for allele in mat pat
  do
ln -s ../${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_${allele}reads_map2ref_DistanceTomaxTSN.bed  .
  done
done

# examine if there is a error
for Head in HT KD SK
do 
cat ${Head}_matReads_patReads_TSS_maxTSNs_ratio0.5-2_map5_*reads.bed | awk '{OFS="\t"} ($6 !=$12) {print $0}'
done


# merge mat, pat, ide into a table
for Head in HT KD SK
do 
intersectBed -s -wao -a <(intersectBed -s -wao -a ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_matreads_map2ref_DistanceTomaxTSN.bed -b  ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_patreads_map2ref_DistanceTomaxTSN.bed \
| cut -f 1-7,14) -b ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_IDEreads_map2ref_DistanceTomaxTSN.bed  | awk '{OFS="\t"} ($16 !=0 ){print $1,$2,$3,$4,$5,$6,$7,$8, $15} ($16 ==0 ){print $1,$2,$3,$4,$5,$6,$7,$8, "NA"}'\
> ${Head}_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed
done

#ln -s /workdir/sc2457/F1_Tissues/Pause_SingleBaseRunOn/pause_Indel/P.CAST_M.B6_F1hybrid.snps.bed .
ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/P.CAST.EiJ_M.C57BL.6J_*aternal_all.fa* .

j=Tissues3_EarlyPause_1bpapart_KSfdr0.1
d=10
d=30
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_mat_patSeq.bed 
 cat ${j}_+-${d}_mat_patSeq.bed  | awk '{OFS="\t"} (substr($4,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8, $9} 
 (substr($4,1,1)=="P") {print  $1,$2,$3,$4,$5, $6, $7, $9, $8}' > ${j}_+-${d}_Early_LateAlleleSeq.bed 

j=Tissues3_LatePause_1bpapart_KSfdr0.1
d=10
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_mat_patSeq.bed 
 cat ${j}_+-${d}_mat_patSeq.bed  | awk '{OFS="\t"} (substr($4,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8, $9} 
 (substr($4,1,1)=="P") {print  $1,$2,$3,$4,$5, $6, $7, $9, $8}' > ${j}_+-${d}_Early_LateAlleleSeq.bed 


j=Tissues3_EarlyPause_BG
d=10
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_mat_patSeq.bed 
 cat ${j}_+-${d}_mat_patSeq.bed  | awk '{OFS="\t"} (substr($4,1,1)=="M") {print $1,$2,$3,$4,$5, $6, $7, $8, $9} 
 (substr($4,1,1)=="P") {print  $1,$2,$3,$4,$5, $6, $7, $9, $8}' > ${j}_+-${d}_Early_LateAlleleSeq.bed 

j=combine_maxPause_noduplicate
d=30
 bedtools getfasta -s -fi mm10.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t"} {print $1, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_mm10.txt 
 paste ${j}.bed  ${j}_mm10.txt > ${j}_+-${d}_mm10_Seq.bed 
 

j=combine_maxPause_noduplicate
d=50
 bedtools getfasta -s -fi mm10.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t"} {print $1, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_mm10.txt 
 paste ${j}.bed  ${j}_mm10.txt > ${j}_+-${d}_mm10_Seq.bed 


# examine the SNPs in the 42 sub_df_SNP.bed
ln -s /workdir/sc2457/F1_Tissues/Pause_SingleBaseRunOn/pause_Indel/P.CAST_M.B6_F1hybrid.snps.bed
j=sub_df_SNP
d=0
 bedtools getfasta -s -fi mm10.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t"} {print $1, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_mm10.txt 
 paste ${j}.bed  ${j}_mm10.txt > ${j}_+-${d}_mm10_Seq.bed 
intersectBed -wao -a ${j}_+-${d}_mm10_Seq.bed -b P.CAST_M.B6_F1hybrid.snps.bed  #manually examine using excel

