# use star aligner

cd /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/
cat P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf | awk '{OFS="\t"; c="_maternal"} {print $1c,$2,$3,$4,$5,$6,$7,$8,$9,$12}' > P.CAST_M.B6_F1hybrid.snps.forSTAR.mat.vcf

cd /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam
mkdir STAR_index_P.CAST.EiJ_M.C57BL.6J_paternal
mkdir STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal


export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH

#withou annotation
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--genomeFastaFiles *_P.CAST.EiJ_M.C57BL.6J_maternal.fa

####### 
cd /workdir/sc2457/F1_Tissues/RNA-seq/
mkdir STAR_LV
cd STAR_LV
ln -s ../LV*fastq.gz .

FASTQ=LV_MB6
export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH

# map1 without annotation
#mat B6 mm10
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ}_BOTH_RNA.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_mat1 \
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

# STAR three stage mapping
# map1 without annotation
# map2 with sj list from map1, mm10 (mat) + vcf files
# map3 b6(mat) with sj list from mat2v, add wasp

#vW:i:1 means alignment passed WASP filtering, and all other values mean it did not pass:
samtools view LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.bam |grep vW:i:1 -c
#14,380,436

f=LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.bam
samtools view $f |grep  vW:i:1 | grep vA:B:c,1 |wc -l
# CAST 7,264,929
samtools view $f |grep  vW:i:1 | grep vA:B:c,2 |wc -l
# B6   7,093,872
samtools view $f |grep  vW:i:1 | grep vA:B:c,3 |wc -l
# No match either allele 21,635

j=`echo $f|rev|cut -d . -f 2- |rev`
rm -r ${j}
mkdir ${j}
samtools view $f |grep  vW:i:1 | awk '{OFS="\t"}{print $0 >> "'${j}'/"$20 }' &

cd $j
rm cast b6
# get the allele specific reads 
# remove chrX and chrY
for v in `ls vA:B:c,1* |grep -v 2 | grep -v 3`
do cat $v |grep -v X_maternal|grep -v Y_maternal | grep -v  MT_maternal >> cast
done

for v in `ls vA:B:c,2* |grep -v 1 | grep -v 3`
do cat $v |grep -v X_maternal|grep -v Y_maternal | grep -v  MT_maternal  >> b6
done

 wc -l cast b6
   7214263 cast
   7008568 b6
  14222831 total


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

wait
cat header cast | samtools sort -@ 20 - > ../${j}.cast.sorted.bam &
cat header b6 | samtools sort -@ 20 - > ../${j}.b6.sorted.bam &

cd ..
samtools view ${j}.b6.sorted.bam |cut -f 3 |uniq > chrOrder.txt
samtools view -H ${j}.b6.sorted.bam |grep @SQ|sed 's/@SQ\tSN:\|LN://g' > genome.txt
# -s	Force “strandedness”
# -split	Treat “split” BAM (i.e., having an “N” CIGAR operation) or BED12 entries as distinct BED intervals.
# -g	Specify a genome file the defines the expected chromosome order in the input files for use with the -sorted option.
for p in b6 cast
do 
  # the read coverage of each SNP position
bedtools coverage -s -split -counts -a <(cat P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed| awk '{OFS="\t"; s="+"} {print $0, ".", s}') -b ${j}.${p}.sorted.bam > ${j}.${p}.plus.snpReadCounts.bed &
bedtools coverage -s -split -counts -a <(cat P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed| awk '{OFS="\t"; s="-"} {print $0, ".", s}') -b ${j}.${p}.sorted.bam > ${j}.${p}.minus.snpReadCounts.bed &
done

wc -l ${j}.*.snpReadCounts.bed
  17482864 LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.minus.snpReadCounts.bed
  17482864 LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.plus.snpReadCounts.bed
  17482864 LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.minus.snpReadCounts.bed
  17482864 LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.plus.snpReadCounts.bed
  69931456 total

for s in plus minus
do 
echo -e "chrm\tsnppos\tmat_allele_count\tpat_allele_count" > ${j}.AlleleHMM_input_${s}.txt
paste ${j}.b6.${s}.snpReadCounts.bed ${j}.cast.${s}.snpReadCounts.bed | awk '{OFS="\t"} ($7+1 !=1|| $14+1 != 1){print $1, $3, $7, $14}' >> ${j}.AlleleHMM_input_${s}.txt
done

ln -s /workdir/sc2457/F1_Tissues/AlleleHMM/AlleleHMM.py .
#j=LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out
# use alleleHMM to identify allelic bias block from RNA-seq data
python2 AlleleHMM.py -p ${j}.AlleleHMM_input_plus.txt -m ${j}.AlleleHMM_input_minus.txt -o LV_MB6_BOTH_RNA -
# HERE

# switch plus amd minus strand in AlleleHMM blocks
for f in LV_MB6_BOTH_RNA_minus_regions_*.bed
do 
	echo $f
    b=`echo $f|rev| cut -d . -f 2- |rev`
    echo $b
	cat $f | awk '{OFS="\t"} ($4!= "S") {print $1, $2, $3, $4, $5, "+"}' > ${b}_filtered_plus.bed &
done
for f in LV_MB6_BOTH_RNA_plus_regions_*.bed
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
	# switch plus and minus strand
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
# switch plus and minus strand
bedtools genomecov -split -bg -ibam ${j}.bam -strand -  >  ${j}.all_plus.bedGraph &
bedtools genomecov -split -bg -ibam ${j}.bam -strand +  > ${j}.all_minus.bedGraph &
wait
k=${j}.all
echo $k

wait
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



ln -s /local/storage/projects/sc2457/F1_Tissues/PolyA_Allele-specific-manuscript/LV_AT_4tunitIntersectNativeHMM* .
wc -l LV_AT_4tunitIntersectNativeHMM_*
  1150 LV_AT_4tunitIntersectNativeHMM_AlleleHMM.bed
  1150 LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed
  1150 LV_AT_4tunitIntersectNativeHMM_tunits.bed

# the intersect between AT window and AlleleHMM blocks
# -s  Force “strandedness”
t=9
intersectBed -u -s -a LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed -b <( cat AlleleHMM_Blocks/LV_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed AlleleHMM_Blocks/LV_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed ) > LV_t1E-0${t}_AT_AlleleHMM.bed #| wc -l 
#intersectBed -v -s -a LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed -b <( cat ../LV_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed ../LV_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed ) | wc -l 
#intersectBed -wao -s -a LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed -b <( cat ../LV_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed ../LV_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed ) 


# AlleleHMM blocks from RNA-seq was just upstream of AT window
# strandness -s, ignore downstream -id
# -fu	Choose first from features in B that are upstream of features in A.
t=9
bedtools closest -s -D a -id -fu -t first -a <(sort-bed LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed) \
-b <( cat AlleleHMM_Blocks/LV_MB6_BOTH_RNA_plus_regions_t1E-0${t}_filtered_minus.bed AlleleHMM_Blocks/LV_MB6_BOTH_RNA_minus_regions_t1E-0${t}_filtered_plus.bed | sort-bed -) \
| awk '{OFS="\t"} ($13+0 < -1 && $13+0 >-100000) {print $0}' \
> LV_t1E-0${t}_AT_AlleleHMM_within10K.bed 
# col 1-6 LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed
# col 7-12 AlleleHMM_Blocks from LV_MB6_BOTH_RNA

cat LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed  | awk '{OFS="\t"} ($6=="+") {print $1, $2, $6} ($6=="-") {print $1, $3, $6}' | sort |uniq |wc -l
#1126
cat LV_t1E-0${t}_AT_AlleleHMM.bed  LV_t1E-0${t}_AT_AlleleHMM_within10K.bed | awk '{OFS="\t"} ($6=="+") {print $1, $2, $6} ($6=="-") {print $1, $3, $6}'  | sort |uniq |wc -l
#462 #          intersect ones                       with AlleleHMM block upstream <10Kb

cat gencode.vM25.annotation.gtf | awk 'BEGIN{OFS="\t"} $3=="gene" {split($10,a,"\""); split($14,b,"\""); print $1, $4-1, $5, a[2], b[2], $7}' > gencode.vM25.annotation.gene.bed
# use the tunits that contain the AT window of interests to identify geneID

intersectBed -wao -s -a LV_AT_4tunitIntersectNativeHMM_tunits.bed -b gencode.vM25.annotation.gene.bed  \
> LV_AT_4tunitIntersectNativeHMM_tunits_gencode.vM25.annotation.geneID.bed
# col 1-6 LV_AT_4tunitIntersectNativeHMM_tunits.bed
# col 7-12 gencode.vM25.annotation.gene.bed
# col 13 intersect length 

cat LV_AT_4tunitIntersectNativeHMM_tunits_gencode.vM25.annotation.geneID.bed | cut -f 10 |sort |uniq > geneID_withATwindow

# -a use the tunit # use the tunits to identify the geneID within that tunits 
# -b LV_t1E-0${t}_AT_AlleleHMM.bed use the LV_AT_4tunitIntersectNativeHMM_intersectRegion.bed that insect with RNA-seq alleleHMM blocks
intersectBed -u -s -a LV_AT_4tunitIntersectNativeHMM_tunits_gencode.vM25.annotation.geneID.bed -b <(cat LV_t1E-0${t}_AT_AlleleHMM.bed LV_t1E-0${t}_AT_AlleleHMM_within10K.bed|cut -f 1-6 | sort | uniq ) | cut -f 10 |sort |uniq \
> geneID_withATwindow.with.nearby.RNA.AlleleHMM.blocks

mv geneID_withATwindow LV_geneID_withATwindow
mv geneID_withATwindow.with.nearby.RNA.AlleleHMM.blocks LV_geneID_withATwindow.with.nearby.RNA.AlleleHMM.blocks

#####
# Use all the AlleleHMM blocks within the Tunits, instead of limiting to 10K upstream of AT window



#####



f=LV_MB6_BOTH_RNA_mat3waspSJ.out.tab


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

# calculate the stability of mRNA
# reads counts of proseq in gene body / read counts of mRNA-seq in exon

# switch the strand for gencode.vM25.annotation.mat.gtf and gencode.vM25.annotation.gtf, 
# because mRNA were mapped to the opposite strand
cat gencode.vM25.annotation.mat.gtf | awk 'BEGIN{FS="\t"; OFS="\t"}($7=="+") {print $1,$2,$3,$4,$5,$6,"-",$8,$9} ($7=="-") {print $1,$2,$3,$4,$5,$6,"+",$8,$9} ($7!="+" && $7 != "-") {print $0}' > gencode.vM25.annotation.mat_strandswitched.gtf &
cat gencode.vM25.annotation.gtf | awk 'BEGIN{FS="\t"; OFS="\t"}($7=="+") {print $1,$2,$3,$4,$5,$6,"-",$8,$9} ($7=="-") {print $1,$2,$3,$4,$5,$6,"+",$8,$9} ($7!="+" && $7 != "-") {print $0}' > gencode.vM25.annotation_strandswitched.gtf &

 ln -s ../STAR_BN/gencode.vM25.*strandswitched.gtf .

# htseq-count [options] alignment_file gff_file
b=LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out
#samtools sort -O SAM -n -@ 20 ${b}.bam > ${b}.nsorted.sam
htseq-count --stranded yes -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation.mat_strandswitched.gtf > ${b}.exon.read.count &
htseq-count --stranded yes -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation.mat_strandswitched.gtf > ${b}.gene.read.count &

b=LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6
#samtools sort -O SAM -n -@ 20 ${b}.sorted.bam > ${b}.nsorted.sam
htseq-count --stranded yes -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation_strandswitched.gtf > ${b}.nsorted.exon.read.count &
htseq-count --stranded yes -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation_strandswitched.gtf > ${b}.nsorted.gene.read.count &


b=LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast
#samtools sort -O SAM -n -@ 20 ${b}.sorted.bam > ${b}.nsorted.sam
htseq-count --stranded yes -t exon -f sam ${b}.nsorted.sam gencode.vM25.annotation_strandswitched.gtf > ${b}.nsorted.exon.read.count &
htseq-count --stranded yes -t gene -f sam ${b}.nsorted.sam gencode.vM25.annotation_strandswitched.gtf > ${b}.nsorted.gene.read.count &


# determine the gene thats in 