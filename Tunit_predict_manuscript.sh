# use 5 state HMM 
for t in BN HT  SK  SP  LG  LV  GI  ST
#for t in HT  SK  SP  LG GI  ST
do
Rscript run.hmm.h5_F1bedgraph.R ${t}_all ${t}_all.dREG.peak.score.bed.gz ${t}_all_plus.bw ${t}_all_minus.bw &
#Rscript run.hmm.h3_F1bedgraph.R ${t}_all ${t}_all.dREG.peak.score.bed.gz ${t}_all_plus.bw ${t}_all_minus.bw &
done

t=KD
Rscript run.hmm.h5_F1bedgraph.R ${t}_all ${t}_all.dREG.peak.score.bed.gz ${t}_all_plus.bw ${t}_all_minus.bw &


t=BN
c=MB6
ci=mat
c=PB6

for t in BN HT  SK  SP  LG  LV  GI  ST
do
for c in MB6 PB6
do
for ci in mat pat
do
#Rscript run.hmm.h3_F1bedgraph.R ${t}_${c}_${ci} ${t}_all.dREG.peak.score.bed.gz ${t}_${c}_all_R1.${ci}_1bp_plus.bw ${t}_${c}_all_R1.${ci}_1bp_minus.bw &
Rscript run.hmm.h5_F1bedgraph.R ${t}_${c}_${ci} ${t}_all.dREG.peak.score.bed.gz ${t}_${c}_all_R1.${ci}_1bp_plus.bw ${t}_${c}_all_R1.${ci}_1bp_minus.bw &
done



# not used. 
# found in human.hg19.h3.log

parallel.sh -j 2 -r "Rscript run.hmm.h3_F1bedgraph.R * BN_all.dREG.peak.score.bed.gz BN_all_plus.bw BN_all_minus.bw chr22.params.Rdata" `Rscript ../../scripts/get.chroms.R H-U_plus.hg19.bw | grep -v "^chr22$" | grep -v "_"`



cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN {OFS="\t"} ($6=="+") { print $0}' > gencode.vM20.annotation_transcript_plus.bed
cat gencode.vM20.annotation_transcript.bed | awk 'BEGIN {OFS="\t"} ($6=="-") { print $0}' > gencode.vM20.annotation_transcript_minus.bed

# select the transcripts in gencode annotation that are likely to expressed
ln -s ~/bin/tuSelecter/* .
# make input annotation file for tuSelecter
./build_txtable.sh gencode.vM20.annotation.gtf > vM20.txtable.out
./tuSelecter.R -p 5 -o gencode_tuSelecter_BN/  vM20.txtable.out BN_all_plus.bw BN_all_minus.bw

