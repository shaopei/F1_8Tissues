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
  # This function take strandness into account
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
  # This function take strandness into account
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
            PREFIX=${bed_dir}/${Head}_${cross}_all_R1   #map2ref_bed/ST_PB6_all_R1
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

# identify imprinted tunits, remove those. do binomial test again with pool reads.
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
    #BinomialTest_IDE only keep block with at lease 1 allele-specific read
done
wait


# bias fdr<=0.1, control fdr > 0.9

for Head in BN LV 
    do 
    	cat ${Head}_TunitProteinSrainEffect_binomtest_fdr.bed | awk 'BEGIN {OFS="\t"} ($11+0 <= 0.1) {print $1,$2,$3,$4,$5, $10, $6,$7,$8,$9, $11}' > ${Head}_TunitProteinSrainEffect_binomtest_fdr0.1.bed
    	cat ${Head}_TunitProteinSrainEffect_binomtest_fdr.bed | awk 'BEGIN {OFS="\t"} ($11+0 > 0.9) {print $1,$2,$3,$4,$5, $10, $6,$7,$8,$9, $11}' > ${Head}_TunitProteinSrainEffect_binomtest_fdr0.9.bed
        cat ${Head}_TunitProteinSrainEffect_binomtest_fdr.bed | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5, $10, $6,$7,$8,$9, $11}' > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed
# chrom, chromStart, chromEnd, HMM, HMM_BinomailTest, strand, mat(B6)Reads, patReads, IdeReads, Binomial pvalue, Binomial fdr

done


## identify the closet genes, including overlap one 
#1.  identify tunits that overlap (in opposite strand)
# -S	Require opposite strandedness
for fdr in fdr0.1 fdr0.9
do 
bedtools closest -S -d -a ${Head}_TunitProteinSrainEffect_binomtest_${fdr}.bed -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
| awk 'BEGIN {OFS="\t"} ($23+0 ==0) {print $0}'  > ${Head}_TunitProteinSrainEffect_binomtest_${fdr}_adjacentTunit.bed
#2. identify the closet one that do not overlap regardless of strandness
                             # exclude those with overlap in the opposite strand to avoid duplicates
# -io	Ignore features in B that overlap A. That is, we want close, yet not touching features only.
bedtools closest -io -d -a <(intersectBed -v -a ${Head}_TunitProteinSrainEffect_binomtest_${fdr}.bed -b ${Head}_TunitProteinSrainEffect_binomtest_${fdr}_adjacentTunit.bed) -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
>> ${Head}_TunitProteinSrainEffect_binomtest_${fdr}_adjacentTunit.bed
done


###Given a gene with AT window,  is the adjacent gene more likely to be biased? ###
# exclude the paired downstream adjancent tunits that are baised to the same direction (i.e. Both B6 or Both CAST) AND on the same strand
cd /workdir/sc2457/F1_Tissues/transcription_level_analysis/AllelicTermination_StainEffectDomain
# from PolyA_Allele-specific-manuscript_2021updates
# tunits that DO overlap with the gene transcript (gencode.vM25.annotation_transcript_protein_coding) with dREG sites 
# exclude chrX and chrY
# Tunits
ln -s /workdir/sc2457/F1_Tissues/PolyA_Allele-specific-manuscript/tunit_protein_coding .
# Tunits that are under strain effect (no evidence of significantly imprinted), with binomial test results
ln -s /workdir/sc2457/F1_Tissues/transcription_level_analysis/domains_cluster_more_than_chance_or_not_tunit_protein/*_TunitProteinSrainEffect_binomtest_fdrAll.bed .
# chrom, chromStart, chromEnd, HMM, HMM_BinomailTest, strand, mat(B6)Reads, patReads, IdeReads, Binomial pvalue, Binomial fdr


# AT windows
ln -s  /workdir/sc2457/F1_Tissues/PolyA_Allele-specific-manuscript/*_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed .
# col 1-6 AT window
# col 7-12 AlleleHMM blocks
# col 13-18 Tunit 
# col 19 B6 or CAST which strain AT bias toward


for Head in BN LV
do
# identify the SrainEffect tunits without AT window
# -b  col 13-18 Tunit
intersectBed -v -s -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
-b <(cat ${Head}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed| cut -f 13-18)| uniq > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow.bed
done

for Head in BN LV
do
# identify the SrainEffect tunits with AT window
# -b  col 13-18 Tunit, col 1-6 AT window
intersectBed -s -wo -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
-b <(cat ${Head}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed| awk 'BEGIN {OFS="\t"} {print $13,$14,$15,$16,$17,$18, $1,$2,$3,$4,$5,$6}' )\
 |awk 'BEGIN {OFS="\t"; a="_AT"} {print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11,$18, $19,$20,$21a,$22,$23}' | uniq \
 > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed
done
# ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed
# col 1-11 ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed
# col 12-17 AT window


## identify the closet genes, including overlap one. if tie, Report the first tie that occurred in the B file.
# exclude the paired downstream adjancent tunits that are baised to the same direction (i.e. Both B6 or Both CAST) AND on the same strand

#1.  identify tunits that overlap (in opposite strand)
# -S	Require opposite strandedness
# SrainEffect tunits with AT window
for Head in BN LV
do
bedtools closest -t first -S -d -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed \
-b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
| awk 'BEGIN {OFS="\t"} ($29+0 ==0) {print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11,$18, $19,$20,$21,$22,$23, $24, $25, $26,$27,$28, $29, $12,$13,$14,$15,$16,$17}' \
 > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunit.bed
# col 1-11 ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed with a pair on opposite strad
# col 12-22 the opposite strad pair in ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed
# col 23 distance
# col 24-29 associated AT window

#2. identify the closet one that do not overlap regardless of strandness
                             # exclude those with overlap in the opposite strand to avoid duplicates
bedtools closest -t first -io -d -a <(intersectBed -v -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunit.bed) \
-b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11,$18, $19,$20,$21,$22,$23,$24, $25, $26,$27,$28, $29, $12,$13,$14,$15,$16,$17}' \
>> ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunit.bed
done


# use paired Tunits only from the opposite strand
for Head in BN LV
do
bedtools closest -t first -S -d -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed \
-b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11,$18, $19,$20,$21,$22,$23, $24, $25, $26,$27,$28, $29, $12,$13,$14,$15,$16,$17}' \
 > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunitOppositeStrand.bed
# col 1-11 ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed with a pair on opposite strad
# col 12-22 the opposite strad pair in ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed
# col 23 distance
# col 24-29 associated AT window
done



#SrainEffect tunits without AT window
## identify the closet genes, including overlap one. if tie, Report the first tie that occurred in the B file.
#1.  identify tunits that overlap (in opposite strand)
# -S	Require opposite strandedness
for Head in BN LV
do
bedtools closest -t first -S -d -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow.bed \
-b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
| awk 'BEGIN {OFS="\t"} ($23+0 ==0) {print $0}'  > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_adjacentTunit.bed
# col 1-11 ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow.bed with a pair on opposite strad
# col 12-22 the opposite strad pair in ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed

#2. identify the closet one that do not overlap regardless of strandness
                             # exclude those with overlap in the opposite strand to avoid duplicates
bedtools closest -t first -io -d -a <(intersectBed -v -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow.bed -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_adjacentTunit.bed) \
-b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
>> ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_adjacentTunit.bed
done



# use paired Tunits only from the opposite strand
for Head in BN LV
do
bedtools closest -t first -S -d -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow.bed \
-b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
 > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_adjacentTunitOppositeStrand.bed
# col 1-11 ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow.bed with a pair on opposite strad
# col 12-22 the opposite strad pair in ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed
done


# compare genes (use protein coding Tunits here) overlap with AT window, and genes overlap with gene(another tunit) without AT window.
# without AT window
# for gene1 without AT window, identify the gene2 with TSS overlap with the 2nd half of the gene1  using R
for Head in BN LV
do
intersectBed -wo -S -a ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow.bed -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
> ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_IntersectTunitOppositeStrand.bed
done

#with AT window
# identify gene2 overlap with AT window
for Head in BN LV
do
cat ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed | awk 'BEGIN {OFS="\t"} {print $12,$13,$14,$15,$16,$17, $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11}' > temp.bed
intersectBed -wo -S -a temp.bed -b ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed \
| awk 'BEGIN {OFS="\t"} {print $7, $8, $9, $10, $11, $12,$13,$14,$15,$16,$17, $18, $19,$20,$21,$22,$23, $24, $25, $26,$27,$28, $29, $1,$2,$3,$4,$5,$6}' \
 > ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunitOppositeStrand.bed
# col 1-11 ${Head}_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow.bed with a pair on opposite strad
# col 12-22 the opposite strad pair in ${Head}_TunitProteinSrainEffect_binomtest_fdrAll.bed
# col 23 distance
# col 24-29 associated AT window
done



# from Find_consistent_blocks_v3.bsh
# strain effect domain
ln -s /workdir/sc2457/F1_Tissues/ImprintingOrGenetics/Combined_MB6andPB6/T8_2Strand_p0.05_effect_strain.bed_cluster .
# seperate strain effect domain to two group, one with AT window, one without

for Head in BN LV
do
# identify the SrainEffect domains without AT window
# -b  col 13-18 Tunit
intersectBed -v -a <(cat T8_2Strand_p0.05_effect_strain.bed_cluster |cut -f 3-| grep ${Head}) \
-b <(cat ${Head}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed| cut -f 13-18| uniq)| uniq > ${Head}_T8_2Strand_p0.05_effect_strain.bed_cluster_withoutATwindow.bed
done

for Head in BN LV
do
# identify the SrainEffect domains with AT window
# -b  col 13-18 Tunit
intersectBed -wa -a <(cat T8_2Strand_p0.05_effect_strain.bed_cluster |cut -f 3-| grep ${Head}) \
-b <(cat ${Head}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed| cut -f 13-18| uniq) | uniq > ${Head}_T8_2Strand_p0.05_effect_strain.bed_cluster_withATwindow.bed
done

 # wc -l *T8_2Strand_p0.05_effect_strain.bed_cluster*
 #   326 BN_T8_2Strand_p0.05_effect_strain.bed_cluster_withATwindow.bed
 #   407 BN_T8_2Strand_p0.05_effect_strain.bed_cluster_withoutATwindow.bed
 #   712 LV_T8_2Strand_p0.05_effect_strain.bed_cluster_withATwindow.bed
 #  1132 LV_T8_2Strand_p0.05_effect_strain.bed_cluster_withoutATwindow.bed
 #  3466 T8_2Strand_p0.05_effect_strain.bed_cluster
 #  6043 total

# use gencode.vM25.annotation.gtf  or Tunits over lap with protein coding gene
ln -s /workdir/sc2457/F1_Tissues/ImprintingOrGenetics/Combined_MB6andPB6/GeneAnnotationInCluster/gencode.vM25.annotation_geneMerged.bed

f=gencode.vM25.annotation_geneMerged.bed

for Head in BN LV
do
	for state in withoutATwindow withATwindow
do
	f=tunit_protein_coding/${Head}_all_h5.preds.full_inProtein_coding.bed
bedtools intersect -wo -a <(cat ${f} |cut -f 1-3) -b  ${Head}_T8_2Strand_p0.05_effect_strain.bed_cluster_${state}.bed |cut -f 4-7 |sort |uniq -c > ${Head}_T8_2Strand_p0.05_effect_strain.bed_cluster_${state}_geneCount.txt
done
done




