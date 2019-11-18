# use 5 state HMM 
for t in BN HT  SK  SP  KD  LV  GI  ST
do
Rscript run.hmm.h5_F1bedgraph.R ${t}_all ${t}_all.dREG.peak.score.bed.gz ${t}_all_plus.bw ${t}_all_minus.bw &
done

# double checked in /workdir/sc2457/F1_Tissues/bigWig/tunit_h5_prediction_replicates