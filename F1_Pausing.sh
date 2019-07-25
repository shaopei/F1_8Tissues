
PREFIX=map2ref_bed/BN_MB6_all_R1
MAT_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
PAT_READ_BED=${PREFIX}.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz

bedtools intersect -a <(zcat BN_all.dREG.peak.score.bed.gz) -b <(zcat ${MAT_READ_BED})

map2ref_bed/BN_MB6_all_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz

for Head in BN HT  SK  SP  LG  LV  GI  ST
do 
echo $Head
bed_dir=map2ref_bed
  for f in ${bed_dir}/${Head}_MB6_*_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz  #each samples from MB6 of the same tissue
    do PREFIX=`echo $f|cut -d . -f 1`
    echo $PREFIX
    MAT_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
    PAT_READ_BED=${PREFIX}.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.sorted.bed.gz
    #IDENTICAL_READ_BED=${PREFIX}.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.sorted.bed.gz

    P=`echo $PREFIX| cut -d / -f 2`
    ln -s HMM_bed/combined_cross/T8_HMM_plus.bed ${P}_HMM_plus.bed
    ln -s HMM_bed/combined_cross/T8_HMM_minus.bed ${P}_HMM_minus.bed
    BinomialTest ${P}_HMM_plus.bed ${P}_HMM_plus ${MAT_READ_BED} ${PAT_READ_BED} &
    BinomialTest ${P}_HMM_minus.bed ${P}_HMM_minus ${MAT_READ_BED} ${PAT_READ_BED} &
    wait
  done

### Ignore above , start from here!!! ####
# identify trasncript annotation that have a dREG sites near the 5' 100bp regions
cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN{OFS="\t"}  ($6=="-") {print $1, $3-100, $3, $4, $5, $6, $0}; 
($6=="+") {print $1, $2, $2+100, $4, $5, $6, $0}' > gencode.vM20.annotation_transcript_100bp.bed

for Head in BN HT  SK  SP  LV  GI  ST KD
do 
intersectBed -a gencode.vM20.annotation_transcript_100bp.bed -b <(zcat Browser/${Head}_all.dREG.peak.score.bed.gz) > ${Head}_gencode.vM20.annotation_transcript_100bp.bed
#|awk 'BEGIN{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14}' 
done

# keep the first 30bp
for Head in BN HT  SK  SP  LV  GI  ST KD
do 
cat ${Head}_gencode.vM20.annotation_transcript_100bp.bed | awk 'BEGIN{OFS="\t"}  ($6=="-") {print $1, $3-30, $3, $4, $5, $6, $0}; 
($6=="+") {print $1, $2, $2+30, $4, $5, $6, $0}'| LC_ALL=C sort -k1,1V -k2,2n --parallel=30 > ${Head}_gencode.vM20.annotation_transcript_30bp.bed 
done


# Keep TRX with SNPs in the first 30bp
unfiltered_snp=/workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered

for Head in BN HT  SK  SP  LV  GI  ST KD
do 
intersectBed -sorted -u -a ${Head}_gencode.vM20.annotation_transcript_30bp.bed -b <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }') > ${Head}_gencode.vM20.annotation_transcript_30bp_wSNP.bed &
#|awk 'BEGIN{OFS="\t"} {print $7,$8,$9,$10,$11,$12,$13,$14}' 
done



echo "remove reads that DONOT overlape with a SNP in *_AMBremoved_sorted_specific.map2ref.sorted.bed"
intersectBed -sorted -u -a <(LC_ALL=C sort -k1,1V -k2,2n --parallel=30 ${MATBOWTIE}_AMBremoved_sorted_specific.map2ref.bed |awk '{print "chr"$0}') \
-b <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }') |gzip > ${MATBOWTIE}_AMBremoved_sorted_specific.map2ref.sorted.bed.gz &
intersectBed -sorted -u -a <(LC_ALL=C sort -k1,1V -k2,2n --parallel=30 ${PATBOWTIE}_AMBremoved_sorted_specific.map2ref.bed |awk '{print "chr"$0}') \
-b <(cat ${unfiltered_snp} |awk '{OFS="\t"}{print "chr"$1, $2-1, $2, $6 }') |gzip > ${PATBOWTIE}_AMBremoved_sorted_specific.map2ref.sorted.bed.gz &
cat ${MATBOWTIE}_AMBremoved_sorted_identical.map2ref.bed | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 |awk '{print "chr"$0}' |gzip > ${MATBOWTIE}_AMBremoved_sorted_identical.map2ref.sorted.bed.gz





bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand + > ${TMPDIR}/$j\_plus.bedGraph
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand - > ${TMPDIR}/$j\_minus.noinv.bedGraph
# calculate coverage in the ${Head}_gencode.vM20.annotation_transcript.bed
bedtools coverage -a $f -b ${MAT_READ_BED} -s | awk 'BEGIN {OFS="\t"; t="_"} {print $1t$2t$3, $0}' |LC_ALL=C sort -k1,1V -k2,2n > ${j}.mat_cov.bed &


bedtools intersect -u -a  KD_PB6_F5_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.1bp.sorted.diff.bed -b /workdir/sc2457/F1_Tissues/Pause_SingleBaseRunOn/KD_gencode.vM20.annotation_transcript.bed  > KD_PB6_F5_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_identical.map2ref.1bp.sorted.diff_inTRX.bed