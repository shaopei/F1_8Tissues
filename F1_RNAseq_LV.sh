
####### 
mkdir STAR_LV
cd STAR_LV
ln -s ../LV*fastq.gz .

FASTQ=LV_MB6
export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH

# map1 without annotation
#mat
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ}_BOTH_RNA.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_mat1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes All

# paternal 
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_paternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ}_BOTH_RNA.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_pat1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes All

#map2 with sj list from map1
#mat with vcf
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ}_BOTH_RNA.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_mat2v \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes All \
--sjdbFileChrStartEnd ${FASTQ}_mat1SJ.out.tab \
--varVCFfile /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.snps.forSTAR.mat.vcf


# paternal 
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_paternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ}_BOTH_RNA.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_pat2 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes All \
--sjdbFileChrStartEnd ${FASTQ}_pat1SJ.out.tab 


#mat3 with sj list from mat2v, wasp
FASTQ=LV_MB6 

#for f in ${FASTQ}*.fastq.gz
#do echo $f
#j=`echo $f |rev |cut -d . -f 3-|rev`
#echo $j
#echo \

f=LV_MB6_BOTH_RNA.fastq.gz
j=LV_MB6_BOTH_RNA

STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn $f \
--runThreadN 40 \
--outFileNamePrefix ${j}_mat3wasp \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI AS nM NM MD jM jI MC ch vA vG vW \
--sjdbFileChrStartEnd ${FASTQ}_mat2vSJ.out.tab \
--varVCFfile /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.snps.forSTAR.mat.vcf \
--waspOutputMode SAMtag
#done


samtools view BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.bam |grep vW:i:1 -c
# 11,852,744

f=BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.bam
samtools view $f |grep  vW:i:1 | grep vA:B:c,1 |wc -l
# CAST  2,698,423
samtools view $f |grep  vW:i:1 | grep vA:B:c,2 |wc -l
# B6    2,653,715
samtools view $f |grep  vW:i:1 | grep vA:B:c,3 |wc -l
# No match either allele 8,560

j=`echo $f|rev|cut -d . -f 2- |rev`
rm -r ${j}
mkdir ${j}
samtools view $f |grep  vW:i:1 | awk '{OFS="\t"}{print $0 >> "'${j}'/"$20 }' &

cd $j
rm cast b6
# remove chrX
for v in `ls vA:B:c,1* |grep -v 2 | grep -v 3`
do cat $v |grep -v X_maternal|grep -v Y_maternal | grep -v  MT_maternal >> cast
done

for v in `ls vA:B:c,2* |grep -v 1 | grep -v 3`
do cat $v |grep -v X_maternal|grep -v Y_maternal | grep -v  MT_maternal  >> b6
done

 wc -l cast b6
   5785623 cast
   5763455 b6
  11549078 total

for f in cast b6
do 
sed -i 's/10_maternal/chr10/g' $f &
done

for c in {11..19}
do  
	for f in cast b6
	do echo sed -i \'s/${c}_maternal/chr${c}/g\' $f "&"
	done
	echo "wait"
done

for c in {1..9}
do  
	for f in cast b6
	do echo sed -i \'s/${c}_maternal/chr${c}/g\' $f "&"
	done
	echo "wait"
done

samtools view -H ../${j}.bam > header
for c in X Y MT
do
echo sed -i \'s/${c}_maternal/chr${c}/g\' header
done


cat header cast | samtools sort -@ 20 - > ../${j}.cast.sorted.bam
cat header b6 | samtools sort -@ 20 - > ../${j}.b6.sorted.bam

cd ..
samtools view ${j}.b6.sorted.bam |cut -f 3 |uniq > chrOrder.txt
samtools view -H ${j}.b6.sorted.bam |grep @SQ|sed 's/@SQ\tSN:\|LN://g' > genome.txt
# -s	Force “strandedness”
# -split	Treat “split” BAM (i.e., having an “N” CIGAR operation) or BED12 entries as distinct BED intervals.
# -g	Specify a genome file the defines the expected chromosome order in the input files for use with the -sorted option.
for p in b6 cast
do 
bedtools coverage -s -split -counts -a <(cat P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed| awk '{OFS="\t"; s="+"} {print $0, ".", s}') -b ${j}.${p}.sorted.bam > ${j}.${p}.plus.snpReadCounts.bed &
bedtools coverage -s -split -counts -a <(cat P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed| awk '{OFS="\t"; s="-"} {print $0, ".", s}') -b ${j}.${p}.sorted.bam > ${j}.${p}.minus.snpReadCounts.bed &
done

wc -l ${j}.*.snpReadCounts.bed
  17482864 BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.minus.snpReadCounts.bed
  17482864 BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.plus.snpReadCounts.bed
  17482864 BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.minus.snpReadCounts.bed
  17482864 BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.plus.snpReadCounts.bed

for s in plus minus
do 
echo -e "chrm\tsnppos\tmat_allele_count\tpat_allele_count" > ${j}.AlleleHMM_input_${s}.txt
paste ${j}.b6.${s}.snpReadCounts.bed ${j}.cast.${s}.snpReadCounts.bed | awk '{OFS="\t"} ($7+1 !=1|| $14+1 != 1){print $1, $3, $7, $14}' >> ${j}.AlleleHMM_input_${s}.txt
done

ln -s /workdir/sc2457/F1_Tissues/AlleleHMM/AlleleHMM.py .
#j=BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out
python2 AlleleHMM.py -p ${j}.AlleleHMM_input_plus.txt -m ${j}.AlleleHMM_input_minus.txt -o BN_MB6_BOTH_RNA -

# switch plus amd minus strand
for f in BN_MB6_BOTH_RNA_minus_regions_*.bed
do 
	echo $f
    b=`echo $f|rev| cut -d . -f 2- |rev`
    echo $b
	cat $f | awk '{OFS="\t"} ($4!= "S") {print $1, $2, $3, $4, $5, "+"}' > ${b}_filtered_plus.bed &
done
for f in BN_MB6_BOTH_RNA_plus_regions_*.bed
do 
		echo $f
    b=`echo $f|rev| cut -d . -f 2- |rev`
    echo $b
	cat $f | awk '{OFS="\t"} ($4!= "S") {print $1, $2, $3, $4, $5, "-"}' > ${b}_filtered_minus.bed &
done





# make bigwig files
export mouse_genome=/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz
export mouse_chinfo=/local/storage/data/mm10/mm10.chromInfo

for p in cast b6
do
	# switch plus amd minus strand
bedtools genomecov -split -bg -ibam ${j}.${p}.sorted.bam -strand - |LC_COLLATE=C sort -k1,1 -k2,2n >  ${j}.${p}_plus.bedGraph &
bedtools genomecov -split -bg -ibam ${j}.${p}.sorted.bam -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}.${p}_minus.bedGraph &
done
wait
for p in cast b6
do
k=${j}.${p}
echo $k
  cat $k\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $k\_minus.inv.bedGraph 
  bedGraphToBigWig $k\_minus.inv.bedGraph ${mouse_chinfo} $k\_minus.bw &
  bedGraphToBigWig $k\_plus.bedGraph ${mouse_chinfo} $k\_plus.bw &
done

# all reads
# switch plus amd minus strand
bedtools genomecov -split -bg -ibam ${j}.bam -strand -  >  ${j}.all_plus.bedGraph &
bedtools genomecov -split -bg -ibam ${j}.bam -strand +  > ${j}.all_minus.bedGraph &
wait
k=${j}.all
echo $k

for c in {10..19}
do 
for f in $k\_plus.bedGraph $k\_minus.bedGraph
do
sed -i 's/'${c}'_maternal/chr'${c}'/g' $f &
done
wait
done

for c in {1..9}
do 
for f in $k\_plus.bedGraph $k\_minus.bedGraph
do
sed -i 's/'${c}'_maternal/chr'${c}'/g' $f &
done
wait
done

for c in X Y MT
do 
for f in $k\_plus.bedGraph $k\_minus.bedGraph
do
sed -i 's/'${c}'_maternal/chr'${c}'/g' $f &
done
wait
done


  cat $k\_minus.bedGraph | grep -v chrMT| LC_COLLATE=C sort -k1,1 -k2,2n| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $k\_minus.inv.bedGraph &
  mv $k\_plus.bedGraph $k\_tempplus.bedGraph
  cat $k\_tempplus.bedGraph | grep -v chrMT| LC_COLLATE=C sort -k1,1 -k2,2n > $k\_plus.bedGraph &
  wait
  rm $k\_tempplus.bedGraph
  bedGraphToBigWig $k\_minus.inv.bedGraph ${mouse_chinfo} $k\_minus.bw &
  bedGraphToBigWig $k\_plus.bedGraph ${mouse_chinfo} $k\_plus.bw &


# only chr10

grep 10_maternal ${j}.all_plus.bedGraph > ${j}_chr10_plus.bedGraph
k=${j}_chr10
c=10
f=$k\_plus.bedGraph 
sed -i 's/'${c}'_maternal/chr'${c}'/g' $f 

 mv $k\_plus.bedGraph $k\_tempplus.bedGraph
  cat $k\_tempplus.bedGraph | grep -v chrMT| LC_COLLATE=C sort -k1,1 -k2,2n > $k\_plus.bedGraph 
  wait
  rm $k\_tempplus.bedGraph
  bedGraphToBigWig $k\_plus.bedGraph ${mouse_chinfo} $k\_plus.bw 


# the intersect between AT window and AlleleHMM blocks
# -s	Force “strandedness”

wc -l BN_AT_4tunitIntersectNativeHMM_intersectRegion.bed
for t in {1..9}
do
echo $t
intersectBed -u -s -a BN_AT_4tunitIntersectNativeHMM_intersectRegion.bed -b <( cat ../BN_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed ../BN_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed ) > BN_t1E-0${t}_AT_AlleleHMM.bed #| wc -l 
intersectBed -v -s -a BN_AT_4tunitIntersectNativeHMM_intersectRegion.bed -b <( cat ../BN_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed ../BN_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed ) | wc -l 
done

t=9
intersectBed -wao -s -a BN_AT_4tunitIntersectNativeHMM_intersectRegion.bed -b <( cat ../BN_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed ../BN_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed ) 


# AlleleHMM blocks from RNA-seq was just upstream of AT window
# strandness -s, ignore downstream -id
# -fu	Choose first from features in B that are upstream of features in A.
t=9
bedtools closest -s -D a -id -fu -t first -a <(sort-bed BN_AT_4tunitIntersectNativeHMM_intersectRegion.bed) -b <( cat ../BN_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed ../BN_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed | sort-bed -) \
| awk '{OFS="\t"} ($13+0 < -1 && $13+0 >-100000) {print $0}' > BN_t1E-0${t}_AT_AlleleHMM_within10K.bed

cat BN_AT_4tunitIntersectNativeHMM_intersectRegion.bed  | awk '{OFS="\t"} ($6=="+") {print $1, $2, $6} ($6=="-") {print $1, $3, $6}' | sort |uniq |wc -l
cat BN_t1E-0${t}_AT_AlleleHMM.bed  BN_t1E-0${t}_AT_AlleleHMM_within10K.bed | awk '{OFS="\t"} ($6=="+") {print $1, $2, $6} ($6=="-") {print $1, $3, $6}'  | sort |uniq |wc -l
#          intersect ones                       with AlleleHMM block upstream <10Kb

# identify the tunits that contain the AT window of interests
intersectBed -wao -s -a BN_AT_4tunitIntersectNativeHMM_tunits.bed -b ../gencode.vM25.annotation.gene.bed  > BN_AT_4tunitIntersectNativeHMM_tunits_gencode.vM25.annotation.geneID.bed
cat BN_AT_4tunitIntersectNativeHMM_tunits_gencode.vM25.annotation.geneID.bed | cut -f 10 |sort |uniq > geneID_withATwindow
#intersectBed -s -a ../gencode.vM25.annotation.gene.bed -b BN_AT_4tunitIntersectNativeHMM_tunits.bed > gencode.vM25.annotation.geneWithATwindow.bed
# use the tunits to identify the geneID within that tunits
intersectBed -u -s -a BN_AT_4tunitIntersectNativeHMM_tunits_gencode.vM25.annotation.geneID.bed -b <(cat BN_t1E-0${t}_AT_AlleleHMM.bed BN_t1E-0${t}_AT_AlleleHMM_within10K.bed|cut -f 1-6 | sort | uniq ) | cut -f 10 |sort |uniq \
> geneID_withATwindow.with.nearby.RNA.AlleleHMM.blocks




f=BN_MB6_BOTH_RNA_mat3waspSJ.out.tab


for c in {10..19}
do 
sed -i 's/'${c}'_maternal/chr'${c}'/g' $f 
done

for c in {1..9}
do 
sed -i 's/'${c}'_maternal/chr'${c}'/g' $f 
done

for c in X Y MT
do 
sed -i 's/'${c}'_maternal/chr'${c}'/g' $f 
done


# paternal genome
# CrossMap.py bam  <chain_file>  <input.bam> [output_file] [options]
CrossMap.py bam /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/pat2ref.chain \
BN_MB6_pat2Aligned.sortedByCoord.out.bam BN_MB6_pat2Aligned.sortedByCoord.out.pat2ref.bam


# reads that mapped to paternal genome perfectly (can be readswithout SNPs or indel)
samtools view  BN_MB6_pat2Aligned.sortedByCoord.out.bam | grep MD:Z:85 > BN_MB6_pat2Aligned.sortedByCoord.out.MD:Z:85.sam
samtools view  BN_MB6_pat2Aligned.sortedByCoord.out.bam | grep MD:Z:84 > BN_MB6_pat2Aligned.sortedByCoord.out.MD:Z:84.sam
samtools view -H BN_MB6_pat2Aligned.sortedByCoord.out.bam > BN_MB6_pat2Aligned.sortedByCoord.out.pat.perfact.match.sam
cat BN_MB6_pat2Aligned.sortedByCoord.out.MD:Z:85.sam BN_MB6_pat2Aligned.sortedByCoord.out.MD:Z:84.sam >> BN_MB6_pat2Aligned.sortedByCoord.out.pat.perfact.match.sam

wc -l BN_MB6_pat2Aligned.sortedByCoord*sam
    4,844,292 BN_MB6_pat2Aligned.sortedByCoord.out.MD:Z:84.sam
   48,517,317 BN_MB6_pat2Aligned.sortedByCoord.out.MD:Z:85.sam

samtools view BN_MB6_pat2Aligned.sortedByCoord.out.bam -F 4 | wc -l
71,318,642

# all reads mapped to paternal genome
samtools view -H BN_MB6_pat2Aligned.sortedByCoord.out.bam |grep @SQ|sed 's/@SQ\tSN:\|LN://g'  > cast.chromInfo

j=BN_MB6_pat2Aligned.sortedByCoord.out
samtools sort -@ 30 ${j}.bam > ${j}.sorted.bam
	# switch plus amd minus strand
bedtools genomecov -split -bg -ibam ${j}.sorted.bam -strand - |LC_COLLATE=C sort -k1,1 -k2,2n >  ${j}_plus.bedGraph &
bedtools genomecov -split -bg -ibam ${j}.sorted.bam -strand + |LC_COLLATE=C sort -k1,1 -k2,2n > ${j}_minus.bedGraph &
wait


# CrossMap.py bam  <chain_file>  <input.bam> [output_file] [options]
CrossMap.py bed /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/pat2ref.chain \
$k\_plus.bedGraph $k\_pat2ref_plus &
CrossMap.py bed /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/pat2ref.chain \
$k\_minus.bedGraph $k\_pat2ref_minus &
wait

k=${j}
cat $k\_pat2ref_plus| awk '{OFS="\t"} {print "chr"$0}' | grep -v chrMT |LC_COLLATE=C sort -k1,1 -k2,2n > $k\_pat2ref_plus.bedGraph &
cat $k\_pat2ref_minus| awk '{OFS="\t"} {print "chr"$0}' | grep -v chrMT |LC_COLLATE=C sort -k1,1 -k2,2n  > $k\_pat2ref_minus.bedGraph &


export mouse_chinfo=/local/storage/data/mm10/mm10.chromInfo
k=${j}_pat2ref
echo $k
  cat $k\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' >  $k\_minus.inv.bedGraph 
  bedGraphToBigWig $k\_minus.inv.bedGraph ${mouse_chinfo} $k\_minus.bw &
  bedGraphToBigWig $k\_plus.bedGraph ${mouse_chinfo} $k\_plus.bw &
wait


# calculate the stability of mRNA
# reads counts of proseq in gene body / read counts of mRNA-seq in exon

cat gencode.vM25.annotation.gtf | awk 'BEGIN{OFS="\t"} $3=="gene" {split($10,a,"\""); split($14,b,"\""); print $1, $4-1, $5,a[2], b[2],$7}' > gencode.vM25.annotation.gene.bed

# htseq-count [options] alignment_file gff_file
b=BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out
samtools sort -O SAM -n -@ 20 ${b}.bam > ${b}.nsorted.sam
htseq-count --stranded yes -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation.mat.gtf > ${b}.exon.read.count 
htseq-count --stranded yes -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation.mat.gtf > ${b}.gene.read.count 

htseq-count --stranded yes --minaqual 20 -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation.mat.gtf > ${b}.exon.read.count-Q20 &
htseq-count --stranded yes --minaqual 20 -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation.mat.gtf > ${b}.gene.read.count-Q20 &


samtools sort -O SAM -n -@ 20 BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.sorted.bam > BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.nsorted.sam
b=BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6
htseq-count --stranded yes -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation.gtf > ${b}.nsorted.exon.read.count 
htseq-count --stranded yes -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation.gtf > ${b}.nsorted.gene.read.count 


b=BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast
samtools sort -O SAM -n -@ 20 ${b}.sorted.bam > ${b}.nsorted.sam
htseq-count --stranded yes -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation.gtf > ${b}.nsorted.exon.read.count 
htseq-count --stranded yes --minaqual 20 -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation.gtf > ${b}.nsorted.exon.read.count-Q20 &
htseq-count --stranded yes -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation.gtf > ${b}.nsorted.gene.read.count 
htseq-count --stranded yes --minaqual 20 -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation.gtf > ${b}.nsorted.gene.read.count-Q20&

# determine the gene thats in 