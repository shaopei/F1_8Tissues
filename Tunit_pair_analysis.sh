### summary of this scrip
##annotate each tunit predict with allele-specific
# map reads from each individual sample from the same tissue to the tunit predicts
# identify blocks share the same allele-bias in MB6 group (A,F,G) or PB6 group (B,C,D,E) of each tissue
# use fishers method to combine p-value of the consistent allele-specificity from multiple samples  
# annotate the tunit predicts with Allelic Bias (AB) at MB6 and PB6, 
# consistent bias to Mat or Pat, and p-value <=0.05. tunit predict without Mat or Pat are labeled as Sym



### annotate each transcript with allele-specific reads mapped to tunit predict
ln -s /workdir/sc2457/F1_Tissues/3rd_batch/map2ref/map2ref_bed/ .

mkdir tunit_preds
cd tunit_preds
ln -s /workdir/sc2457/F1_Tissues/bigWig/*_all_h5.preds.ext.bed .
# seperate tunit ptredicts into plus and minus strand
for f in *.preds.ext.bed
do 
PREFIX=`echo $f|rev|cut -d . -f 3- |rev`
echo $PREFIX
grep plus $f > ${PREFIX}_plus.bed
grep minus $f > ${PREFIX}_minus.bed
done

cd ..

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

  bedtools coverage -s -a <(cat $f | LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}') -b ${MAT_READ_BED} -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $1,$2,$3,$7,$8,$9,$10}' > ${j}.mat_cov.bed &
  bedtools coverage -s -a <(cat $f | LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}') -b ${PAT_READ_BED} -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $1,$2,$3,$7,$8,$9,$10}' > ${j}.pat_cov.bed &
  #bedtools coverage -s -a $f -b ${IDENTICAL_READ_BED} -sorted| awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $1,$2,$3,$7,$8,$9,$10}' |LC_ALL=C sort -k1,1V -k2,2n > ${j}.iden_cov.bed &
  wait

  # filter the block and only keep block with at lease 1 allele-specific read
  join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.5 ${j}.mat_cov.bed ${j}.pat_cov.bed | awk 'BEGIN{OFS="\t"} ($5+$6 >0) {print $0}' | \
  awk 'BEGIN{OFS="\t"; t=","} ($5>$6) {print $2, $3, $4, "M"t$5t$6, $5, $6, "-" }  ($5<$6) {print $2, $3, $4, "P"t$5t$6, $5, $6, "-"}  ($5==$6) {print $2, $3, $4, "S"t$5t$6, $5, $6, "-"} '  \
  | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${j}.merged_cov.bed

  mv ${j}.mat_cov.bed ${j}.pat_cov.bed toremove
  # output of BinomialTestFor_merged_cov.bed.py:(hmm+BinomialTest) if p-value <= 0.05, remain what it got from hmm (can ne M,P, or S), otherwise S.
  python ${PL}/BinomialTestFor_merged_cov.bed.py ${j}.merged_cov.bed ${j}_binomtest.bed
  mv ${j}.merged_cov.bed toremove

  python ${PL}/FalsePosFor_merged_cov.bed.py ${j}_binomtest.bed ${FDR_SIMS} ${FDR_CUTOFF} > ${j}_binomtest_FDR${FDR_CUTOFF}.txt &
  #  awk 'NR==1 { print $0 } NR>1 && $9 <= thresh { print $0 }'  thresh=$(awk 'END {print $6}' ${j}_binomtest_FDR${FDR_CUTOFF}.txt) < ${j}_binomtest.bed  > ${j}_interestingHets.bed
}

# map reads from each individual sample from the same tissue and cross to the tunit predicts

for Head in BN HT SK SP LG LV GI ST
    do 
    echo $Head #BN
    bed_dir=map2ref_bed
    for cross in MB6 PB6
        do echo $cross #MB6
          for f in ${bed_dir}/${Head}_${cross}_*_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz  #each samples from MB6 of the same tissue
            do echo $f  #map2ref_bed/ST_PB6_B_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz
            PREFIX=`echo $f|cut -d . -f 1`   #map2ref_bed/ST_PB6_B_R1
            echo $PREFIX
            MAT_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
            PAT_READ_BED=${PREFIX}.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
            #IDENTICAL_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz

            P=`echo $PREFIX| cut -d / -f 2` #ST_PB6_B_R1
            ln -s tunit_preds/${Head}_all_h5.preds_plus.bed ${P}_plus.bed
            ln -s tunit_preds/${Head}_all_h5.preds_minus.bed ${P}_minus.bed

            if [[ "$cross" == "MB6" ]] ; then
                BinomialTest ${P}_plus.bed ${P}_plus ${MAT_READ_BED} ${PAT_READ_BED} &
                BinomialTest ${P}_minus.bed ${P}_minus ${MAT_READ_BED} ${PAT_READ_BED} &
            elif [[ "$cross" == "PB6" ]]  ; then
            # switch -m and -p for PB6
                BinomialTest ${P}_plus.bed ${P}_plus ${PAT_READ_BED} ${MAT_READ_BED} &
                BinomialTest ${P}_minus.bed ${P}_minus ${PAT_READ_BED} ${MAT_READ_BED} &
            fi
          done
    done
    wait
done

# identify blocks share the same allele-bias in MB6 group (A,F,G) or PB6 group (B,C,D,E)
## use fishers method to combine p-value from multiple samples  
## combine p-value in one file
for t in BN HT SK SP LG LV GI ST
  do 
R --vanilla --slave --args $(pwd) $t < Find_consistent_blocks_v2_tunit.R
done

wc -l *_ABconsistent* |paste - -
sed  -i '1i #chrm\tchrmStart\tchrmEnd\twinP_count_all_samples\tp_value_all_samples\twinP\tp_value_Fisher\tp_value_individual_samples' *_ABconsistent*
# use sed to insert the first line





join_AB(){
  file1=$1
  file2=$2
  file3=$3
  cat ${file1} ${file2} |grep -v "#" | awk 'BEGIN{t="_"} {print $1t$2t$3}' |sort -k1,1|uniq > ${file1}.tmp1
  join -a 1 -t $'\t' -e S -j 1 -o 1.1,2.5 ${file1}.tmp1  <(<${file1} awk 'BEGIN{OFS="\t"; t="_"}  {print $1t$2t$3,$0}' | sort -k1,1) > ${file1}.tmp2
  join -a 1 -t $'\t' -e S -j 1 -o 1.1,1.2,2.5 ${file1}.tmp2 <(<${file2} awk 'BEGIN{OFS="\t"; t="_"} {print $1t$2t$3,$0}' | sort -k1,1) | awk 'BEGIN{OFS="\t"} {split($1,a,"_"); print a[1],a[2],a[3], $2, $3}' |sort-bed - > ${file3}
  rm ${file1}.tmp*
}


# annotate the tunit predicts with Allelic Bias (AB) at MB6 and PB6 
for t in BN HT SK SP LG LV GI ST
    do
for cross in MB6 PB6
do
  for s in plus minus
    do 
    join_in1=tunit_preds/${t}_all_h5.preds_${s}.bed 
    join_in2=${t}_${cross}_tunit_preds_${s}_ABconsistent_FisherMethodP0.05.bed
    join_out=tunit_preds/${t}_all_h5.preds_${s}_${cross}_ABconsistent_FisherMethodP0.05.bed
    join_AB ${join_in1} ${join_in2} ${join_out}
  done
done
done

cd tunit_preds/
# form clusters
cluster_distance=1000000
d=1M
for t in BN HT SK SP LG LV GI ST
    do
    for cross in MB6 PB6
        do 
        #chrm   chrmStart       chrmEnd winP_count_all_samples  p_value_all_samples     winP    p_value_Fisher  p_value_individual_samples
        cat ${t}_all_h5.preds_plus_${cross}_ABconsistent_FisherMethodP0.05.bed |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4, substr($5,1,1), "+"}' |grep -v "#" > ${t}_all_h5.preds_2strands_${cross}_ABconsistent_FisherMethodP0.05.bed
        cat ${t}_all_h5.preds_minus_${cross}_ABconsistent_FisherMethodP0.05.bed |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4, substr($5,1,1),"-"}'|grep -v "#" >> ${t}_all_h5.preds_2strands_${cross}_ABconsistent_FisherMethodP0.05.bed
        bedtools cluster -d ${cluster_distance} -i <(sort-bed ${t}_all_h5.preds_2strands_${cross}_ABconsistent_FisherMethodP0.05.bed) > ${t}_all_h5.preds_2strands_${cross}_ABconsistent_FisherMethodP0.05_cluster${d}.bed
    done
done

# pairs within clusters
# label Con/Dis/OneS and TSS distance
for t in BN HT SK SP LG LV GI ST
    do
    for cross in MB6 PB6
        do 
        python Pair_of_bedregions_withincluster.py ${t}_all_h5.preds_2strands_${cross}_ABconsistent_FisherMethodP0.05_cluster${d}.bed ${t}_${cross}_paires_within_cluster${d}.txt
    done
done

#make histogram
for t in BN HT SK SP LG LV GI ST
    do
    for cross in MB6 PB6
        do 
        Rscript Tunit_pair_analysis_histgram.R  $(pwd) ${t}_${cross}_paires_within_cluster${d}.txt ${t}_${cross}_paires_within_cluster${d}.pdf 
    done
done


#####
# find the tunit that runover the TSS of the other tunit of the pair
for t in BN HT SK SP LG LV GI ST
    do
    for cross in MB6 PB6
        do 
        echo ${t}_${cross} > TSS_run_over_${t}_${cross}_${d}.log 
        python TSS_run_over.py ${t}_${cross}_paires_within_cluster${d}.txt ${t}_${cross}_paires_within_cluster${d}_pairsWithin1M.txt >> TSS_run_over_${t}_${cross}_${d}.log &
    done
  wait
done

maxD=3000
#make histogram
for t in BN HT SK SP LG LV GI ST
    do
    for cross in MB6 PB6
        do 
        Rscript TSS_run_over.R  $(pwd) ${t}_${cross}_paires_within_cluster${d}_pairsWithin1M.txt ${maxD} ${t}_${cross}_paires_within_cluster${d}_pairsWithin1M_${maxD}.pdf ${t}_${cross}_paires_within_cluster${d}_pairsWithin1M_${maxD}RunOver.pdf 
    done
done



#####
# see the fraction of Con and Dis of TSS within AS Tunits
d=1M
for t in BN HT SK SP LG LV GI ST
    do
    for cross in MB6 PB6
        do 
        #echo ${t}_${cross} > TSS_within_Tunit_${t}_${cross}_${d}.log 
        python TSS_within_Tunit.py ${t}_${cross}_paires_within_cluster${d}.txt TSS_within_Tunit_${t}_${cross}_${d} > TSS_within_Tunit_${t}_${cross}_${d}.log &
    done 
done
wait

#rm TSS_within_Tunit_andRunOver_*.log
cat TSS_within_Tunit*log > TSS_within_Tunit_counts_report.txt
Rscript TSS_within_Tunit.R

d=1M
for t in BN HT SK SP LG LV GI ST
    do
    for cross in MB6 PB6
        do 
        #echo ${t}_${cross} > TSS_within_Tunit_${t}_${cross}_${d}.log 
        python TSS_within_Tunit_andRunOver.py ${t}_${cross}_paires_within_cluster${d}.txt TSS_within_Tunit_${t}_${cross}_${d} > TSS_within_Tunit_andRunOver_${t}_${cross}_${d}.log &
    done 
done

cat TSS_within_Tunit_andRunOver_*.log  > TSS_within_Tunit_andRunOver_report.txt
R --vanilla --slave --args $(pwd) TSS_within_Tunit_andRunOver_report.txt TSS_within_Tunit_andRunOver_report.pdf <TSS_within_Tunit_andRunOver.R







#### Below, past code,  don't use


# use intersect in annotate tunit_predict
# get a stats
rm  tunit_aHMM_hit.summary.txt
for t in BN HT SK SP LG LV GI ST
    do 
    for strand in plus minus
        do
        tunit_bed=${t}_all_h5.preds_${strand}.bed
        ahmm_bed=combined_cross/${t}_HMM_${strand}.bed
        log=${tunit_bed}.stats.txt

        #total tunit
        echo -e "${tunit_bed}\t " > ${log}
        wc -l ${tunit_bed} >> ${log}
        # tunit without no interesct --> no allele-specificifty
        bedtools intersect -a ${tunit_bed} -b ${ahmm_bed}  -wao |awk 'BEGIN{OFS="\t"} ($NF==0) {print $0}' | wc -l > tmp
        echo "0" >> tmp
        cat tmp |paste - - >> ${log}
        # tunit with intersect
        bedtools intersect -a ${tunit_bed} -b ${ahmm_bed}  -wao |awk 'BEGIN{OFS="\t"} ($NF!=0) {print $0}' | cut -f 4 |uniq -c  |sort | awk '{print $1}' |uniq -c | awk 'BEGIN{OFS="\t"} {print $1, $2}' >> ${log}
        cat ${log} |awk 'BEGIN{OFS="\t"}  NR<=5 {print $1}'|paste - - - - - >> tunit_aHMM_hit.summary.txt
    done
done

# annotate
for t in BN HT SK SP LG LV GI ST
    do 
    for strand in plus minus
        do
        tunit_bed=${t}_all_h5.preds_${strand}.bed
        ahmm_bed=combined_cross/${t}_HMM_${strand}.bed
        log=${tunit_bed}.stats.txt

        # tunit without no interesct --> no allele-specificifty
        bedtools intersect -a ${tunit_bed} -b ${ahmm_bed}  -wao |awk 'BEGIN{OFS="\t"} ($NF==0) {print $1,$2,$3,$4,$5,$6,$7,$8,$9, "S"}' > ${tunit_bed}_S
        # tunit with 1 intersect
        bedtools intersect -a ${tunit_bed} -b ${ahmm_bed}  -wao |awk 'BEGIN{OFS="\t"} ($NF!=0) {print $0}' | cut -f 4 |uniq -c  |sort | awk '{print $1}' |uniq -c | awk 'BEGIN{OFS="\t"} {print $1, $2}' >> ${log}
        cat ${log} |awk 'BEGIN{OFS="\t"}  NR<=5 {print $1}'|paste - - - - - >> tunit_aHMM_hit.summary.txt


# stats
# How many ASE in each tunit?


# identify tunit pairs


