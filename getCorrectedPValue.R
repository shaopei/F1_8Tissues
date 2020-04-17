#R --vanilla --slave --args $(pwd) input_fp target_col cut_off output < getCorrectedPValue.R
#R --vanilla --slave --args $(pwd) HT_allReads_TSS_maxTSNs_SNPs20bp_binomtest.bed 9 0.1 HT_allReads_TSS_maxTSNs_SNPs20bp_binomtest.fdr0.1R.bed< getCorrectedPValue.R

# get corrected p value from a col of the input file

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_fp=args[2]
target_col=as.numeric(args[3])
cut_off=as.numeric(args[4])
output = args[5]

df <- read.table(input_fp, sep="\t")
colnames(df)=c("chrm","chrmStart", "chrmEnd", "MSP", "MSP_BinomialTest", "mat_allele_count", "pat_allele_count", "identical_reads_count", "Binom_p_value", "strand")
df$fdr = p.adjust(df[,target_col], method = "fdr")

write.table(df[df$fdr<=cut_off,], file=output, quote =F, sep = "\t", row.names = F)
