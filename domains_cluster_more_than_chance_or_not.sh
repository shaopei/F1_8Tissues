# Condition on two groups of TUs: 
# One that has significant allelic specificity. One that does NOT. 
# Ask: Do TUs that are adjacent to a significant change tend to share the same change? 


## two group of TU
cd /workdir/sc2457/F1_Tissues/transcription_level_analysis/domains_cluster_more_than_chance_or_not_tunit_protein
# use tunits overlap with protein
#ln -s /workdir/sc2457/F1_Tissues/bigWig/tunit .  
ln -s /workdir/sc2457/F1_Tissues/PolyA_Allele-specific-manuscript/tunit_protein_coding .

ln -s  /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5 .

### calculate allele specifcity within each tunit predicts
PL=/workdir/sc2457/alleleDB/alleledb_pipeline_mouse
FDR_SIMS=10
FDR_CUTOFF=0.1

BinomialTest(){
  # BinomialTest $f $j $MAT_READ_BED $PAT_READ_BED
  # Need to seperate plus and minus strand before using this function. This function ignore strandness in the output
  f=$1  #bed file containing region to perform test
  j=$2  # prefix for output
  MAT_READ_BED=$3
  PAT_READ_BED=$4

  bedtools coverage -s -a <(cat $f | LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}') -b <(zcat ${MAT_READ_BED} | grep -v chrX |grep -v chrY | grep -v chrMT) -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3t$6, $1,$2,$3,$7,$8,$9,$10, $6}' > ${j}.mat_cov.bed &
  bedtools coverage -s -a <(cat $f | LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}') -b <(zcat ${PAT_READ_BED} | grep -v chrX |grep -v chrY | grep -v chrMT) -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3t$6, $1,$2,$3,$7,$8,$9,$10, $6}' > ${j}.pat_cov.bed &
  #bedtools coverage -s -a $f -b ${IDENTICAL_READ_BED} -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $1,$2,$3,$7,$8,$9,$10}' |LC_ALL=C sort -k1,1V -k2,2n > ${j}.iden_cov.bed &
  wait

  # filter the block and only keep block with at lease 1 allele-specific read
  join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.5,1.9 ${j}.mat_cov.bed ${j}.pat_cov.bed | awk 'BEGIN{OFS="\t"} ($5+$6 >0) {print $0}' | \
  awk 'BEGIN{OFS="\t"; t=","} ($5>$6) {print $2, $3, $4, "M"t$5t$6, $5, $6, "-", $7 }  ($5<$6) {print $2, $3, $4, "P"t$5t$6, $5, $6, "-", $7 }  ($5==$6) {print $2, $3, $4, "S"t$5t$6, $5, $6, "-", $7} '  \
  | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${j}.merged_cov.bed

  mv ${j}.mat_cov.bed ${j}.pat_cov.bed toremove
  # output of BinomialTestFor_merged_cov.bed.py:(hmm+BinomialTest) if p-value <= 0.05, remain what it got from hmm (can ne M,P, or S), otherwise S.
  python2 ${PL}/BinomialTestFor_merged_cov.bed.py ${j}.merged_cov.bed ${j}_binomtest.bed
  mv ${j}.merged_cov.bed toremove
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed ${j}_binomtest_fdr.bed < getFDR.R
  # python2 ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt &
  #  awk 'NR==1 { print $0 } NR>1 && $9 <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_interestingHets.bed
}

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

  # filter the block and only keep block with at lease 1 allele-specific read
  join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.5,1.9 ${j}.mat_cov.bed ${j}.pat_cov.bed > ${j}.merged_cov.temp
  join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.5,1.7 ${j}.merged_cov.temp ${j}.ide_cov.bed \
  | awk 'BEGIN{OFS="\t"} ($5+$6 >0) {print $0}' | \
  awk 'BEGIN{OFS="\t"; t=","} ($5>$6) {print $2, $3, $4, "M"t$5t$6, $5, $6, $7, $8 }  ($5<$6) {print $2, $3, $4, "P"t$5t$6, $5, $6, $7, $8 }  ($5==$6) {print $2, $3, $4, "S"t$5t$6, $5, $6, $7, $8 } '  \
  | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${j}.merged_cov.bed

  mv ${j}.mat_cov.bed ${j}.pat_cov.bed toremove
  # output of BinomialTestFor_merged_cov.bed.py:(hmm+BinomialTest) if p-value <= 0.05, remain what it got from hmm (can ne M,P, or S), otherwise S.
  python2 ${PL}/BinomialTestFor_merged_cov.bed.py ${j}.merged_cov.bed ${j}_binomtest.bed
  mv ${j}.merged_cov.bed toremove
  R --vanilla --slave --args $(pwd) ${j}_binomtest.bed ${j}_binomtest_fdr.bed < getFDR.R
  # python2 ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt &
  #  awk 'NR==1 { print $0 } NR>1 && $9 <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_interestingHets.bed
}

# map reads from each individual sample from the same tissue and cross to the tunit predicts

for Head in BN LV #HT SK SP LV GI ST KD
    do 
    echo $Head #BN
    bed_dir=map2ref_1bpbed_map5
    for cross in MB6 PB6
        do echo $cross #MB6
          #for f in ${bed_dir}/${Head}_${cross}_all_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz
           # do echo $f  #map2ref_bed/ST_PB6_B_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz
            PREFIX=${bed_dir}/${Head}_${cross}_all_R1   #map2ref_bed/ST_PB6_B_R1
            echo $PREFIX
            MAT_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz
            PAT_READ_BED=${PREFIX}.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz
            #IDENTICAL_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz
            P=${Head}_${cross}_all_h5
            ln -s tunit_protein_coding/${Head}_all_h5.preds.full_inProtein_coding.bed ${P}.bed

            if [[ "$cross" == "MB6" ]] ; then
                BinomialTest ${P}.bed ${P} ${MAT_READ_BED} ${PAT_READ_BED} &
            elif [[ "$cross" == "PB6" ]]  ; then
            # switch -m and -p for PB6
                BinomialTest ${P}.bed ${P} ${PAT_READ_BED} ${MAT_READ_BED} &
            fi
    done
done
    wait

# identify imprinted tunits, remove those. do binomial test again witj pool reads.
for Head in BN LV 
    do 
intersectBed -wb -a ${Head}_MB6_all_h5_binomtest_fdr.bed -b ${Head}_PB6_all_h5_binomtest_fdr.bed | \
awk 'BEGIN {OFS="\t"} ($2==$13 && $3==$14 && substr($4,1,1) == substr($15,1,1) && $11+0 <=0.1 && $22+0 <= 0.1 ){print $1,$2,$3,$4,$5,$10}' \
> ${Head}_Tunit_imprinted_effect_fdr0.1.bed
done

for Head in BN LV 
    do intersectBed -v -a tunit_protein_coding/${Head}_all_h5.preds.full_inProtein_coding.bed \
    -b ${Head}_Tunit_imprinted_effect_fdr0.1.bed > ${Head}_all_h5.preds.full_inProtein_coding_NoImprinted.bed
done

# group reads from MB6 and PB6
bed_dir=map2ref_1bpbed_map5
for Head in BN LV #HT SK SP LV GI ST KD
    do 
    zcat ${bed_dir}/${Head}_*_all_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz \
    |grep -v chrX |grep -v chrY | grep -v chrMT| LC_ALL=C sort -k1,1V -k2,2n --parallel=30 |gzip > ${Head}_MAT_READ_BED.temp.gz &
    zcat ${bed_dir}/${Head}_*_all_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted.bed.gz \
    |grep -v chrX |grep -v chrY | grep -v chrMT| LC_ALL=C sort -k1,1V -k2,2n --parallel=30 |gzip > ${Head}_PAT_READ_BED.temp.gz &
    zcat ${bed_dir}/${Head}_*_all_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.map5.1bp.sorted.bed.gz \
    |grep -v chrX |grep -v chrY | grep -v chrMT| LC_ALL=C sort -k1,1V -k2,2n --parallel=30 |gzip > ${Head}_IDE_READ_BED.temp.gz &

done


for Head in BN LV #HT SK SP LV GI ST KD
    do 
    echo $Head #BN
    P=${Head}_TunitProteinSrainEffect
    ln -s ${Head}_all_h5.preds.full_inProtein_coding_NoImprinted.bed ${P}.bed
    BinomialTest_IDE ${P}.bed ${P} ${Head}_MAT_READ_BED.temp.gz ${Head}_PAT_READ_BED.temp.gz ${Head}_IDE_READ_BED.temp.gz &
done
wait


# bias fdr<=0.1, control fdr > 0.9

for Head in BN LV 
    do 
    	cat ${Head}_TunitProteinSrainEffect_binomtest_fdr.bed | awk 'BEGIN {OFS="\t"} ($11+0 <= 0.1) {print $1,$2,$3,$4,$5, $10, $6,$7,$8,$9, $11}' > ${Head}_TunitProteinSrainEffect_binomtest_fdr0.1.bed
    	cat ${Head}_TunitProteinSrainEffect_binomtest_fdr.bed | awk 'BEGIN {OFS="\t"} ($11+0 > 0.9) {print $1,$2,$3,$4,$5, $10, $6,$7,$8,$9, $11}' > ${Head}_TunitProteinSrainEffect_binomtest_fdr0.9.bed
        cat ${Head}_TunitProteinSrainEffect_binomtest_fdr.bed | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5, $10, $6,$7,$8,$9, $11}' > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed
done


#1.  identify tunits that overlap (in opposite strand)
# -S	Require opposite strandedness
for fdr in fdr0.1 fdr0.9
do 
bedtools closest -S -d -a ${Head}_TunitProteinSrainEffect_binomtest_${fdr}.bed -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
| awk 'BEGIN {OFS="\t"} ($23+0 ==0) {print $0}'  > ${Head}_TunitProteinSrainEffect_binomtest_${fdr}_adjacentTunit.bed
#2. identify the closet one that do not overlap
bedtools closest -io -d -a <(intersectBed -v -a ${Head}_TunitProteinSrainEffect_binomtest_${fdr}.bed -b ${Head}_TunitProteinSrainEffect_binomtest_${fdr}_adjacentTunit.bed) -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
>> ${Head}_TunitProteinSrainEffect_binomtest_${fdr}_adjacentTunit.bed
done



############ below old, NOT used ########

# bias fdr<=0.1, control fdr > 0.9
# biased, strain effect, tunits
for Head in BN LV 
    do 
    for strand in plus minus
        do
intersectBed -wb -a ${Head}_MB6_all_h5_${strand}_binomtest_fdr.bed -b ${Head}_PB6_all_h5_${strand}_binomtest_fdr.bed | \
awk 'BEGIN {OFS="\t"} ($2==$12 && $3==$13 && substr($4,1,1) != substr($14,1,1) && $10+0 <=0.1 && $20+0 <= 0.1 ){print $0}' \
> ${Head}_${strand}_Tunit_strain_effect_fdr0.1.temp

intersectBed -wb -a <(cat ${Head}_all_h5_${strand}.bed | cut -f 1-6) \
-b <(cat ${Head}_${strand}_Tunit_strain_effect_fdr0.1.temp | awk 'BEGIN {OFS="\t"} (substr($4,1,1)=="M") {print $1, $2, $3, "B6"} (substr($4,1,1)=="P") {print $1, $2, $3, "CAST"}') \
|awk 'BEGIN {OFS="\t"} ($2==$8 && $3==$9){print $1,$2,$3,$4, $10, $6}' \
> ${Head}_${strand}_Tunit_strain_effect_fdr0.1.bed

# control, UNbiased, tunits
intersectBed -wb -a ${Head}_MB6_all_h5_${strand}_binomtest_fdr.bed -b ${Head}_PB6_all_h5_${strand}_binomtest_fdr.bed | \
awk 'BEGIN {OFS="\t"} ($2==$12 && $3==$13 && $5=="S" && $15=="S" && $10+0 >0.9 && $20+0 >0.9){print $0}' \
> ${Head}_${strand}_Tunit_strain_effect_fdr0.9.temp

intersectBed -wb -a <(cat ${Head}_all_h5_${strand}.bed | cut -f 1-6) \
-b <(cat ${Head}_${strand}_Tunit_strain_effect_fdr0.9.temp | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "Sym"} ') \
|awk 'BEGIN {OFS="\t"} ($2==$8 && $3==$9){print $1,$2,$3,$4, $10, $6}' \
> ${Head}_${strand}_Tunit_strain_effect_fdr0.9.bed

cat ${Head}_${strand}_Tunit_strain_effect_fdr0.1.bed ${Head}_${strand}_Tunit_strain_effect_fdr0.9.bed |sort-bed - > ${Head}_${strand}_Tunit_strain_effect_fdr0.1AND0.9.bed

intersectBed -wb -a ${Head}_MB6_all_h5_${strand}_binomtest_fdr.bed -b ${Head}_PB6_all_h5_${strand}_binomtest_fdr.bed | \
awk 'BEGIN {OFS="\t"} ($2==$12 && $3==$13) {print $0}' \
| awk 'BEGIN {OFS="\t"} ($5 !="S" && substr($4,1,1) =="M" && substr($14,1,1)=="P" && $10+0 <=0.1 && $20+0 <= 0.1 ){print $1, $2, $3, "B6", $10, $20} 
                        ($5 !="S" && substr($4,1,1) =="P" && substr($14,1,1)=="M" && $10+0 <=0.1 && $20+0 <= 0.1 ){print $1, $2, $3, "CAST", $10, $20} 
                        ($5 =="S"){print $1, $2, $3, "Sym", $10, $20}' \
> ${Head}_${strand}_Tunit_strain_effect_fdrALL.temp

intersectBed -wb -a <(cat ${Head}_all_h5_${strand}.bed | cut -f 1-6) -b ${Head}_${strand}_Tunit_strain_effect_fdrALL.temp \
| awk 'BEGIN {OFS="\t"} ($2==$8 && $3==$9){print $1,$2,$3,$4, $10, $6, $11, $12}' > ${Head}_${strand}_Tunit_strain_effect_fdrALL.bed


done
done




