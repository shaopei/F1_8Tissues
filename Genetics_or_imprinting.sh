# genetics, both B6 or both Cast
# imprinting, both M or both P

ln -s /workdir/sc2457/F1_Tissues/3rd_batch/map2ref/*strict* .


for f in *_agreeCount_strict_agreenAmong_*_samples.bed
do echo $f #BN_MB6_all_R1_HMM_minus_agreeCount_strict_agreenAmong_3_samples.bed
head=`echo $f |cut -d _ -f 1-8`
mv $f ${head}.bed
done

# intersect the blocks from two cross MB6 and PB6
for f in *MB6_all_R1_HMM_*_agreeCount_strict.bed
do head=`echo $f |cut -d _ -f 1`
tail=`echo $f |cut -d _ -f 3-`
s=`echo $f |cut -d _ -f 6`

bedtools intersect -wo -f 0.1 -F 0.1 -sorted -a <( cat $f| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4'}) -b <(cat ${head}_PB6_${tail} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}') \
> ${head}_MB6_PB6_${s}_intersect_0.1.bed
bedtools intersect -wo -sorted -a <( cat $f| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4'}) -b <(cat ${head}_PB6_${tail} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}') \
> ${head}_MB6_PB6_${s}_intersect_1bp.bed

# bedtools intersect -wo -sorted -a <( cat $f| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4'}) -b <(cat ${head}_PB6_${tail} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}') \
# |awk 'BEGIN{OFS="\t"} ($NF >10 ){print $0}' > ${head}_MB6_PB6_${s}_intersect_10.bed
done


for f in *_intersect_*.bed
do 
cat $f | awk 'BEGIN {OFS="\t"} substr($4,1,1)==substr($8,1,1) {print $0}' > ${f}_strain.bed # M.P in PB6 is reversed
cat $f | awk 'BEGIN {OFS="\t"} substr($4,1,1)!=substr($8,1,1) {print $0}' > ${f}_imprint.bed
done


diff LV_MB6_PB6_minus_intersect_0.1.bed_imprint.bed LV_MB6_PB6_minus_intersect_1bp.bed_imprint.bed | awk '{print $2":"$3"-"$4,$6":"$7"-"$8, $0}'^C
[sc2457@cbsudanko ImprintingOrGenetics]$ head LV_MB6_PB6_minus_intersect_0.1.bed_imprint.bed | awk '{print $1":"$2"-"$3,$5":"$6"-"$7, $0}'


ln -s /workdir/sc2457/F1_Tissues/3rd_batch/map2ref/*allbut1* .
for f in *_agreeCount_allbut1_agreenAmong_*_samples.bed
do echo $f #BN_MB6_all_R1_HMM_minus_agreeCount_strict_agreenAmong_3_samples.bed
head=`echo $f |cut -d _ -f 1-8`
mv $f ${head}.bed
done

# intersect the blocks from two cross MB6 and PB6
for f in *MB6_all_R1_HMM_*_agreeCount_allbut1.bed
do head=`echo $f |cut -d _ -f 1`
tail=`echo $f |cut -d _ -f 3-`
s=`echo $f |cut -d _ -f 6`

bedtools intersect -wo -f 0.1 -F 0.1 -sorted -a <( cat $f| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4'}) -b <(cat ${head}_PB6_${tail} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}') \
> ${head}_MB6_PB6_${s}_intersect_0.1.bed
bedtools intersect -wo -sorted -a <( cat $f| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4'}) -b <(cat ${head}_PB6_${tail} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}') \
> ${head}_MB6_PB6_${s}_intersect_1bp.bed
