library(BioVenn)
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/TID_TSR_maxTSS")
Head="BN"
a=read.table(paste(Head,"_temp_A", sep=""))
b=read.table(paste(Head,"_temp_B", sep=""))
c=read.table(paste(Head,"_temp_C", sep=""))
biovenn <- draw.venn(a$V1, b$V1, c$V1, subtitle=Head, nrtype="abs")
#https://cran.r-project.org/web/packages/BioVenn/vignettes/BioVenn.html

library(eulerr)
#set.seed(1)
#BN
combo <- c(A = 3225, B = 1093, C = 2324, "A&B" = 476, "A&C" = 410, "B&C" = 126, "A&B&C"=94)
plot(euler(combo), main="BN")
#LV
combo <- c(A = 7227, B = 1482, C = 10113, "A&B" = 745, "A&C" = 2133, "B&C" = 313, "A&B&C"=242)
plot(euler(combo), main="LV")
