tfbs.drawLogo(tfs, file.pdf="logos.pdf", motif_id=c("M00216_2.00", "M09433_2.00", "M00216_2.00") )

file.twoBit <- "/home/sc2457/mm10.2bit"

## Example 1: Scan the whole genome
## Scan 2bit file within whole genome to find motif binding site
r1.scan <- tfbs.scanTFsite( tfs, file.twoBit, ncores = 20);


write.table(r1.scan$result[[1]], file="mm10-TBP_M00216_2.00.bed",
            quote=F, sep="\t", row.names = F,
            col.names = T)

write.table(r1.scan$result[[2]], file="mm10-TBP_M09433_2.00.bed",
            quote=F, sep="\t", row.names = F,
            col.names = T)



########
## Attach the package firstly
library(rtfbsdb);
## Create db from downloaded dataset
db<- CisBP.download("Mus_musculus");

## Query the CisBP dataset and select the motifs for a transcription factor of interest
tfs<- tfbs.createFromCisBP(db, family_name="TBP");
file.twoBit <- "/home/sc2457/mm10.2bit"
#r1.scan <- tfbs.scanTFsite( tfs, file.twoBit, ncores = 20, threshold = 3);

Head="LV" #"BN"
b6.TATA.Bed = read.table(paste(Head,"_allReads_TSS_maxTSNs_binomtest_-35to-20_mat.bed",sep=""))
CAST.TATA.Bed = read.table(paste(Head,"_allReads_TSS_maxTSNs_binomtest_-35to-20_pat.bed",sep=""))


file.B6.fa = "P.CAST.EiJ_M.C57BL.6J_maternal_all.fa"
file.CAST.fa = "P.CAST.EiJ_M.C57BL.6J_paternal_all.fa"

b6.scan <- tfbs.scanTFsite( tfs, file.B6.fa, b6.TATA.Bed, ncores = 20, threshold = 3);
dim(b6.TATA.Bed)
dim(b6.scan$result[[1]])

b6.6.scan <- tfbs.scanTFsite( tfs, file.B6.fa, b6.TATA.Bed, ncores = 20, threshold = 6);

CAST.scan <- tfbs.scanTFsite( tfs, file.CAST.fa, CAST.TATA.Bed, ncores = 20, threshold = 3);


b6.M09433 = b6.scan$result[[2]]
temp <- data.frame(do.call(rbind, strsplit(as.character(b6.M09433$chrom), "_")))
b6.M09433  = cbind(b6.M09433,temp[1])

b6.M00216 = b6.scan$result[[1]]
temp <- data.frame(do.call(rbind, strsplit(as.character(b6.M00216$chrom), "_")))
b6.M00216 = cbind(b6.M00216,temp[1])

CAST.M09433 = CAST.scan$result[[2]]
temp <- data.frame(do.call(rbind, strsplit(as.character(CAST.M09433$chrom), "_")))
CAST.M09433 = cbind(CAST.M09433,temp[1])

CAST.M00216 = CAST.scan$result[[1]]
temp <- data.frame(do.call(rbind, strsplit(as.character(CAST.M00216$chrom), "_")))
CAST.M00216 = cbind(CAST.M00216,temp[1])


df.M09433 = merge(b6.M09433, CAST.M09433, 
	by=c("chromStart","chromEnd", "X1","name", "strand"),
	all = FALSE,
	suffixes = c(".b6",".cast"),
	sort = F,
	)
dim(df.M09433)
dim(b6.M09433)
dim(CAST.M09433)

df.M09433$score.ave = (df.M09433$score.b6 + df.M09433$score.cast)/2

df.M00216 = merge(b6.M00216, CAST.M00216, 
	by=c("chromStart","chromEnd", "X1","name", "strand"),
	all = FALSE,
	suffixes = c(".b6",".cast"),
	sort = F,
	)
dim(df.M00216)
dim(b6.M00216)
dim(CAST.M00216)
df.M00216$score.ave = (df.M00216$score.b6 + df.M00216$score.cast)/2




write.table(df.M09433[,c(3,1,2,7,11,5,14)], file=paste(Head,"_allReads_TSS_maxTSNs_binomtest_motifM09433.bed", sep = ""),
            quote=F, sep="\t", row.names = F,
            col.names = T)

write.table(df.M00216[,c(3,1,2,7,11,5,14)], file=paste(Head,"_allReads_TSS_maxTSNs_binomtest_motifM00216.bed", sep = ""),
            quote=F, sep="\t", row.names = F,
            col.names = T)
