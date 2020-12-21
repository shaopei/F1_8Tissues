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
# CAST 3,444,200  2,803,362
samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam | grep vA:B:c,2 |wc -l
# B6 2,592,491  3,019,146
samtools view BN_MB6_F23_RNA_test32waspAligned.sortedByCoord.out.bam |grep vA:B:c,3 |wc -l
# No match either allele 376,516


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



