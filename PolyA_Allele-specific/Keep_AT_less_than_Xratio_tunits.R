#R --vanilla --slave --args $(pwd) ${Head}_AT_3tunitIntersectNativeHMM_intersectRegion.bed ${Head}_AT_3tunitIntersectNativeHMM_tunits.bed ${Head}_AT_4tunitIntersectNativeHMM_intersectRegion.bed ${Head}_AT_4tunitIntersectNativeHMM_tunits.bed 0.5 < Keep_AT_less_than_Xratio_tunits.R

args=(commandArgs(TRUE))

folder = args[1]
f1_p = args[2] #"BN_AT_3tunitIntersectNativeHMM_intersectRegion.bed"
f2_p = args[3] #"BN_AT_3tunitIntersectNativeHMM_tunits.bed"
f3_p = args[4]
f4_p = args[5]
r = as.numeric(args[6])

setwd(folder)

tunit <-  read.table(f2_p, header = F)
AT <- read.table(f1_p, header = F)
AT.tunit.length.ratio <- (AT$V3-AT$V2)/(tunit$V3-tunit$V2)


write.table(AT[AT.tunit.length.ratio  <= as.numeric(r),], file = f3_p, quote = FALSE, sep = "\t",
	row.names = FALSE, col.names = FALSE)

write.table(tunit[AT.tunit.length.ratio  <= as.numeric(r),], file = f4_p, quote = FALSE, sep = "\t",
	row.names = FALSE, col.names = FALSE)