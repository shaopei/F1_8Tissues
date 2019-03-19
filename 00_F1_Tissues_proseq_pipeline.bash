### Organize the raw sequences files
# Run proseqHT_multiple_adapters_sequencial.bsh to get QC reads

for f in /local/storage/projects/RawSequenceFiles/2019-01-04-F1-mouse/raw_data/*_1.fq.gz
do j=`echo $f |rev |cut -d _ -f 2- |rev`
jj=`echo $j |rev | cut -d '/' -f 1 |rev`
ln -s ${j}_1.fq.gz ${jj}_R1.fastq.gz
ln -s ${j}_2.fq.gz ${jj}_R2.fastq.gz
done

export mouse_genome=/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz
export mouse_chinfo=/local/storage/data/mm10/mm10.chromInfo

BWAIDX=/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz
CHINFO=/local/storage/data/mm10/mm10.chromInfo

output=proseqHT_0
mkdir ${output}
bash proseqHT_multiple_adapters_sequencial.bsh -I \*.fastq.gz -i $mouse_genome -c $mouse_chinfo -T  ${output} -O  ${output} 


mkdir In_sample_list
while read p; do
  echo "$p"
  mv ${p}* In_sample_list/.
done <sample_list.txt

cd In_sample_list
for f in *.sort.bam
do
	j=`echo ${f} |rev|cut -d . -f 3- |rev `
	echo $j
samtools sort ${f} -@ 30 -o ${j}.Csort.bam
java -jar /programs/bin/picard-tools/picard.jar CollectInsertSizeMetrics \
      I=${j}.Csort.bam \
      O=${j}_insert_size_metrics.txt \
      H=${j}_insert_size_histogram.pdf 
done

echo 'sample_ID' > insert_size.txt
cat PoolB_Lung_D_Liver_AGTCA_insert_size_metrics.txt | awk 'NR==7 {print $0}' >>insert_size.txt
for f in *.sort.bam
do
	j=`echo ${f} |rev|cut -d . -f 3- |rev `
	echo ${j} >> insert_size.txt
	cat ${j}_insert_size_metrics.txt | awk 'NR==8 {print $0}' >>insert_size.txt
done

cat insert_size.txt | paste - - > insert_size_results.txt


while read field1 field2 ; do
	echo mv ${field1}_insert_size_histogram.pdf ${field2}_${field1}_insert_size_histogram.pdf
done <sample_group.txt


### seperate the reads before duplicates, so that I can estimate the duplicating rate of each sample
for name in `ls ${TMPDIR}/noadapt/*.fastq | awk -F"/" '{print $NF}' | rev | cut -d \. -f 2- | cut -d _ -f2- | rev| sort | uniq`
   do
   #gzip -d ${TMPDIR}/noadapt/${name}_R1.fastq.gz  &
   #gzip -d ${TMPDIR}/noadapt/${name}_R2.fastq.gz  &
   #wait
   ## Separate into distinct fastq files.  Also trims off specified lengths of sequence.
   # use the sepIndex.py without Trims 3-prime UMI and J-barcode 
   python sepIndex_beforeDeduplicate.py ${TMPDIR}/noadapt/${name}_R1.fastq ${TMPDIR}/noadapt/${name}_R2.fastq ${TMPDIR}/sep_withdups/${name} | tee ${TMPDIR}/${name}_sep_withdups.log &
 done

wait
for f in ${TMPDIR}/noadapt/*.fastq
do gzip $f &
done

for f in *noadapt_sep_withdups.log
do 
name=`echo $f | rev |cut -d _ -f 4- |rev`
cat $f | awk -v name="${name}" '{OFS="_"} NR> 6 {print name, $0}' | cut -d " " -f 1-2 |  sed "s/://" >> sep_withdedup_read_counts.txt
done

cat sep_withdedup_read_counts.txt | sort > sep_withdedup_read_counts.sort
cat ../sample_list.txt |sort > sample_list.sort
join -1 1 -2 1 -o 1.1,2.2 sample_list.sort sep_withdedup_read_counts.sort


## seperate the raw reads nased on J-barcode
for name in `ls *.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 2- | cut -d _ -f2- | rev| sort | uniq`
   do
   python sepIndex_beforeDeduplicate.py ${name}_R1.fastq.gz ${name}_R2.fastq.gz sep_raw/${name} | tee sep_raw/${name}_sep_raw.log &
 done

rm sep_withdedup_read_counts.txt
for f in *_sep_raw.log
do 
name=`echo $f | rev |cut -d _ -f 3- |rev`
cat $f | awk -v name="${name}" '{OFS="_"} NR> 6 {print name, $0}' | cut -d " " -f 1-2 |  sed "s/://" >> sep_withdedup_read_counts.txt
done

cat sep_withdedup_read_counts.txt | LC_ALL=C sort > sep_withdedup_read_counts.sort
cat sample_list.txt |LC_ALL=C sort > sample_list.sort
LC_ALL=C join -1 1 -2 1 -o 1.1,2.2 sample_list.sort sep_withdedup_read_counts.sort


cat rename_sample_list.txt | awk '{OFS="\t"; t="_"; s="-";c="chr"} {print $1, $2t$3}' > rename_sample_pair.txt
cat sample_list.txt  | awk '{OFS="\t"; t1="_*_";t="_"; s="-";c="chr"} {print $2t1$3, $2t$4t$3}' > rename_sample_pair.txt

while read c1 c2 
do
  echo mv ${c1}_R1.fastq.gz ${c2}_R1.fastq.gz >> rename_sample_pair_2.bsh
  echo mv ${c1}_R2.fastq.gz ${c2}_R2.fastq.gz >> rename_sample_pair_2.bsh
done < rename_sample_pair.txt



### QC in the samples
# mouse mm10_GRCm38
# make bed from gtf from GENCODE
gtf2bed < gencode.vM20.annotation.gtf > gencode.vM20.annotation.bed

#transcript
cat gencode.vM20.annotation.bed | cut -f 1-8 |awk 'BEGIN {OFS="\t"} ($8=="transcript"){print $0}' |LC_ALL=C sort --temporary-directory=/workdir/sc2457/tmp/ --parallel=10 -V  > gencode.vM20.annotation_transcript.bed

R --vanilla --slave  < getCounts.R


