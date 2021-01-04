dir_in_sc2457=/workdir/sc2457/F1_Tissues/RNA-seq
cd ${dir_in_sc2457}
ln -s ../alleledb_strandSpecific_9_specifyFASTQ_mouse_SingleEndBowtie.sh .

fastq_dir=/workdir/sc2457/F1_Tissues/RNA-seq
for file in ${fastq_dir}/*_RNA.fastq.gz  #*/
do 
  name=`echo ${file} | rev | cut -d / -f 1 |cut -d . -f 3 | rev`
  echo $name
  cd ${dir_in_sc2457}
echo "bash alleledb_strandSpecific_9_specifyFASTQ_mouse_SingleEndBowtie.sh \
${name} \
PersonalGenome_P.CAST_M.B6 \
/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam \
/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.bed \
${fastq_dir} \
${name} \
/workdir/sc2457/alleleDB/alleledb_pipeline_mouse \
${dir_in_sc2457}/PIPELINE_StrandSpecific_P.CAST_M.B6.mk \
0.1 \
ase \
0 > allelicbias-PersonalGenome_P.CAST_M.B6-${name}.log 2>&1  & " >${name}.bsh
done



# use star aligner

cd /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam

mkdir STAR_index_P.CAST.EiJ_M.C57BL.6J_paternal
mkdir STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal
mkdir STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal_withgencode.vM25.annotation.gtf
mkdir STAR_index_P.CAST.EiJ_M.C57BL.6J_diploid
mkdir STAR_index_GRCm38_68.fa

# make annotation for maternal genome (mm10)
sed -i 's/chr10/10_maternal/g' gencode.vM25.annotation.mat.gtf
for c in {11..19}
do echo sed -i \'s/chr${c}/${c}_maternal/g\' gencode.vM25.annotation.mat.gtf
done
for c in {1..9}
do echo sed -i \'s/chr${c}/${c}_maternal/g\' gencode.vM25.annotation.mat.gtf
done

sed -i 's/chrX/X_maternal/g' gencode.vM25.annotation.mat.gtf
sed -i 's/chrY/Y_maternal/g' gencode.vM25.annotation.mat.gtf
sed -i 's/chrMT/MT_maternal/g' gencode.vM25.annotation.mat.gtf
sed -i 's/MY_maternal/MT_maternal/g' gencode.vM25.annotation.mat.gtf


export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir STAR_index_GRCm38_68.fa \
--genomeFastaFiles GRCm38_68.fa

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal_withgencode.vM25.annotation.gtf \
--genomeFastaFiles *_P.CAST.EiJ_M.C57BL.6J_maternal.fa \
--sjdbGTFfile gencode.vM25.annotation.mat.gtf \
--sjdbOverhang 84

#withou annotation
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--genomeFastaFiles *_P.CAST.EiJ_M.C57BL.6J_maternal.fa

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir STAR_index_P.CAST.EiJ_M.C57BL.6J_paternal \
--genomeFastaFiles *_P.CAST.EiJ_M.C57BL.6J_paternal.fa

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir STAR_index_P.CAST.EiJ_M.C57BL.6J_diploid \
--genomeFastaFiles *_P.CAST.EiJ_M.C57BL.6J_*aternal.fa

FASTQ_PATH=.
FASTQ=BN_MB6_F23_RNA
PL=/workdir/sc2457/alleleDB/alleledb_pipeline_mouse
PGENOME_PATH=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam

zcat ${FASTQ_PATH}/${FASTQ}.fastq.gz | gzip -c | ${PL}/alleledb_filter_input.sh ${PL} - | cat - > ${FASTQ}_temp.fastq
#total 29824806 ns 20189
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal_withgencode.vM25.annotation.gtf \
--readFilesIn ${FASTQ}_temp.fastq \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test2 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard




# try STAR  without using VCF personal variants.
# with annnotation
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal_withgencode.vM25.annotation.gtf \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard

# try STAR  without using VCF personal variants.
# withOUT annnotation
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test3 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard

# paternal 
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_paternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes All

# try STAR  without using VCF personal variants.
# withOUT annnotation BUT add sj.tab from previous run
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test32 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--sjdbFileChrStartEnd ${FASTQ}_test3SJ.out.tab


# paternal 
STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_paternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test42 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes All \
--sjdbFileChrStartEnd ${FASTQ}_test4SJ.out.tab


# try STAR  with using VCF personal variants.
# withOUT annnotation BUT add sj.tab from previous run
cat P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf | awk '{OFS="\t"; c="_maternal"} {print $1c,$2,$3,$4,$5,$6,$7,$8,$9,$12}' > P.CAST_M.B6_F1hybrid.snps.forSTAR.mat.vcf


STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test32v \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI AS nM NM MD jM jI MC ch vA vG \
--sjdbFileChrStartEnd ${FASTQ}_test3SJ.out.tab \
--varVCFfile /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.snps.forSTAR.mat.vcf

STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_maternal \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_test32wasp \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI AS nM NM MD jM jI MC ch vA vG vW \
--sjdbFileChrStartEnd ${FASTQ}_test3SJ.out.tab \
--varVCFfile /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.snps.forSTAR.mat.vcf \
--waspOutputMode SAMtag


samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam |grep  vW:i:1 | grep vA:B:c,1 |wc -l
# CAST  2,698,423
samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam |grep  vW:i:1 | grep vA:B:c,2 |wc -l
# B6    2,653,715
samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam |grep  vW:i:1 | grep vA:B:c,3 |wc -l
# No match either allele 8,560

samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam | grep vA:B:c,1 |wc -l
# CAST 2,803,362
samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam | grep vA:B:c,2 |wc -l
# B6 3,019,146
samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam |grep vA:B:c,3 |wc -l
# No match either allele 10,112


[sc2457@cbsudanko STAR]$ samtools view  BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam |grep vW:i:1 |head
NB500947:967:HCTKLBGXF:1:11306:21307:15427      16      10_maternal     3385350 255     85M     *       0       0       TGTTTATTTATTTATTTATTTATTTATTTGTTTGTTTGTTTGTTTTTATTTTGTTTTTCGAGACATGGTTTTTCTGTGTAGCCCT    /EEEEEAEE/EEEEEEEEAA/AEEAEEEAEE/E<EEEEEEEEEEA/EEEEEEEEEEEE/EEEEEAEEEEEEEEEEEEAEEA//AA   NH:i:1  HI:i:1  AS:i:83 nM:i:0   NM:i:0  MD:Z:85 jM:B:c,-1       jI:B:i,-1       vA:B:c,2        vG:B:i,3385374  vW:i:1
NB500947:967:HCTKLBGXF:1:12307:12340:6857       16      10_maternal     3390108 255     85M     *       0       0       TAGTTGGCTAACTTAATACTGGGGTATCCATTGTCACAGGTTCAGGGAACAATTAGCATCGTGAGCCATTGTGCTGATTCTCCCA    EEEEEEAEEEEEEEEEEEEEEE6EE6EEE6EEEEAEEAEAAAEEEEEE<EEEEEEEEE/EEEEEAEEEEEE/AEEAAEEEA6AA6   NH:i:1  HI:i:1  AS:i:81 nM:i:1   NM:i:1  MD:Z:57G27      jM:B:c,-1       jI:B:i,-1       vA:B:c,1        vG:B:i,3390164  vW:i:1
NB500947:967:HCTKLBGXF:2:23303:3992:16102       16      10_maternal     3390108 255     85M     *       0       0       TAGTTGGCTAACTTAATACTGGGGTATCCATTGTCACAGGTTCAGGGAACAATTAGCATCGTGAGCCATTGTGCTGATTCTCCCA    EEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA   NH:i:1  HI:i:1  AS:i:81 nM:i:1   NM:i:1  MD:Z:57G27      jM:B:c,-1       jI:B:i,-1       vA:B:c,1        vG:B:i,3390164  vW:i:1
NB500947:967:HCTKLBGXF:3:21609:7259:12981       16      10_maternal     3390108 255     85M     *       0       0       TAGTTGGCTAACTTAATACTGGGGTATCCATTGTCACAGGTTCAGGGAACAATTAGCATCGTGAGCCATTGTGCTGATTCTCCCA    EEEE//EEEEEE/EAEEEEEEAAEEEEEEEEEEEEEAE<EEEEEEEAEEEE<EEEEEEEEEA/EAEEEEEEEEEEEEEEEAAA/A   NH:i:1  HI:i:1  AS:i:81 nM:i:1   NM:i:1  MD:Z:57G27      jM:B:c,-1       jI:B:i,-1       vA:B:c,1        vG:B:i,3390164  vW:i:1
NB500947:967:HCTKLBGXF:2:13311:25736:17045      16      10_maternal     3419899 255     85M     *       0       0       TAGGTAGTCCACGAATATGCTTGTTTACAGTAACTTTGGGGGGTTGCCATTAGAATGAAAGCATGAAGCAAAAATTGTTTTTGTA    /E6AAEE/AEAEE<EE/AE/EEEEEEEEE<EAEEEEEAEEEEEAEEEAEEEEEEEEA/EAEEAEEEEAAEEE/EEEEEEEAAAAA   NH:i:1  HI:i:1  AS:i:81 nM:i:1   NM:i:1  MD:Z:46T38      jM:B:c,-1       jI:B:i,-1       vA:B:c,1        vG:B:i,3419944  vW:i:1
NB500947:967:HCTKLBGXF:1:23109:8065:13989       16      10_maternal     3433982 255     85M     *       0       0       TCAGTTAAGGAATCTGCTACGTAGAGAGGCTGCCAAAGTTCCCAAAGCTTGTCAGTGAACAGGTCAGAGTATGCATGGAAACCCA    E/EEEEE<AAAEEEE6EEA<AEEEEE<EEEEEE/<EEEE/EEEAEEEEEEEEEEAE/EAEEEEEEEEEEEEEEEEEEEEEAAAA/   NH:i:1  HI:i:1  AS:i:83 nM:i:0   NM:i:0  MD:Z:85 jM:B:c,-1       jI:B:i,-1       vA:B:c,2        vG:B:i,3434001  vW:i:1
NB500947:967:HCTKLBGXF:4:22401:18603:5464       16      10_maternal     3435714 255     85M     *       0       0       TGAGACTTCGTCTAAAGAGCCATCAGTCTGGGTCTGCAAAAGGAACAGATCTTGACTAAGTAAAGGACTCCACACCCCACAAAAG    EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA   NH:i:1  HI:i:1  AS:i:83 nM:i:0   NM:i:0  MD:Z:85 jM:B:c,-1       jI:B:i,-1       vA:B:c,2        vG:B:i,3435746  vW:i:1

MD:Z:57G27 
samtools view BN_MB6_F23_RNA_test32vAligned.sortedByCoord.out.bam |grep NB500947:967:HCTKLBGXF:1:12307:12340:6857 
samtools view BN_MB6_F23_RNA_test4Aligned.sortedByCoord.out.bam |grep NB500947:967:HCTKLBGXF:1:12307:12340:6857 


## using diploid genome, no annotation
cd /workdir/sc2457/F1_Tissues/RNA-seq/STAR
FASTQ_PATH=.
FASTQ=BN_MB6_F23_RNA
PL=/workdir/sc2457/alleleDB/alleledb_pipeline_mouse
PGENOME_PATH=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam

STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_diploid \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_diploid_test3 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard


STAR --genomeDir /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/STAR_index_P.CAST.EiJ_M.C57BL.6J_diploid \
--readFilesCommand zcat \
--readFilesIn ${FASTQ_PATH}/${FASTQ}.fastq.gz \
--runThreadN 20 \
--outFileNamePrefix ${FASTQ}_diploid_test32 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes All \
--sjdbFileChrStartEnd ${FASTQ}_diploid_test3SJ.out.tab



####### 
mkdir STAR_BN
cd STAR_BN
ln -s ../BN*fastq.gz .

FASTQ=BN_MB6

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
FASTQ=BN_MB6 

for f in ${FASTQ}*.fastq.gz
do echo $f
j=`echo $f |rev |cut -d . -f 3-|rev`
echo $j
#echo \
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
done


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