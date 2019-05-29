


# merge 

### annotate each transcript with consistent AlleleHMM blocks from each tissues from Find_consistent_blocks_v2.bsh
ln -s /workdir/sc2457/F1_Tissues/Find_consistent_blocks/*_ABconsistent_FisherMethodP0.05.bed .

# merge blocks that overlap
for f in *_ABconsistent_FisherMethodP0.05.bed  #BN_MB6_combined_R1_HMM_minus_ABconsistent_FisherMethodP0.05.bed
	do PREFIX=`echo $f|cut -d _ -f 1-2`  #BN_MB6
	s=`echo $f | cut -d _ -f 6` #minus
	bedtools merge -i  $f -c 6 -o distinct -d 0 > ${PREFIX}_${s}_ABconsistent_FMp0.05_merged.bed
done

# seperate tunit ptrficts into plus and minus strand
for f in *.preds.ext.bed
do 
PREFIX=`echo $f|rev|cut -d . -f 3- |rev`
echo $PREFIX
grep plus $f > ${PREFIX}_plus.bed
grep minus $f > ${PREFIX}_minus.bed
done

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


