#R --vanilla --slave --args $(pwd) map5_bw/ HT _allReads_TSS < getMaxTSN_cbsudanko.R
# only use bed regions with at least 5 VALID mat AND 5 pat VALID reads

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
map5_dir=args[2]
Tissue=args[3]
name_body=args[4]

t=Tissue
#map5_dir="map5_bw/"
tss_dir="./"
#name_body="_allReads_TSS"

library("bigWig")
read_read_mat_S <-function (file.plus.bw, file.minus.bw , bed6, step=2, navg = 20, times=1, use.log=FALSE)
{
  bw.plus  <- load.bigWig( file.plus.bw )
  bw.minus <- load.bigWig( file.minus.bw )
  
  hCountMatrix <- bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6[,c(1:6)] , step=step, abs.value=TRUE, op = "sum")
  hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(bed6[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
  if (!use.log){
    hmat <- times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) ;
  } else {
    hmat <- log(times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) + 1) ;
  }
  #avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
  
  unload.bigWig(bw.plus);
  unload.bigWig(bw.minus);
  
  return(hmat);    
}


# all map5 reads
map5.file.plus.bw <- paste(map5_dir,t,"_PB6_F5N6_dedup_QC_end_map5_plus.bw", sep="")
map5.file.minus.bw <- paste(map5_dir,t,"_PB6_F5N6_dedup_QC_end_map5_minus.bw", sep="")

TSS <- read.table(paste(tss_dir,t,name_body,".bed", sep =""), header = F)
names(TSS)[4]="TSNCount"
names(TSS)[5]="ReadsCount"


AT <- TSS[,1:6] 



step=1
#newAT<- NULL
for (i in 1:NROW(AT)){
  # count the reads at each base
  map5.all <- read_read_mat_S (map5.file.plus.bw, map5.file.minus.bw, AT[i,], step) 
  m <- max(map5.all)
  if (m > 0){
    # identify the base with of max map5 reads (maxTSN)
    map5.peaks <- which(map5.all == m)
    # make new table (newAT), the bed region (TSS here) will be duplicated, depends on how many maxTSN in the TSS
    # only TSS with >0 read counts will stay
    for (j in 1:sum(map5.all == m)){
      a = cbind(AT[i,], map5.peaks=map5.peaks[j], maxReadCount = m )
      write.table(a, file=paste(t,name_body,"_maxTSNs",".bed", sep = ""), append = T, quote=F, row.names=F, col.names=F, sep="\t")
      #newAT <- rbind(newAT, a)
    }
  }
}