#R --vanilla --slave --args $(pwd) input_fp target_col cut_off < getHistFromCol.R
#R --vanilla --slave --args $(pwd) HT_allReads_TSS_maxTSNs_DistToDownstreamSNP.bed 11 100 < getHistFromCol.R
# get histpgram from a col of the input file

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_fp=args[2]
target_col=as.numeric(args[3])
cut_off=as.numeric(args[4])

df <- read.table(input_fp, sep="\t")

h <-hist(log10(df[,target_col]+1), 
         breaks = seq(0, ceiling(max(log10(df[,target_col]+1))), 1),
         plot=FALSE)
h$counts=h$counts/sum(h$counts)

pdf(paste(input_fp, "_Full_hist.pdf", sep = ""))
plot(h,col="red" 
     ,density=25     
     ,ylab="Proportion"
     #, xlim=c(0,7)
     ,xlab="Distance to closet downstream SNP (log10)"     
     ,las=2
     #,xaxt='n'
     ,main= input_fp)
dev.off()

pdf(paste(input_fp, "_0-100_hist.pdf", sep = ""))
#l=df[,target_col][df[,target_col]<=100]
#l=l[l>=0]

h2 <-hist(df[,target_col][df[,target_col]<=cut_off],
          breaks = seq(-10,cut_off,1),
          xlim = c(0,cut_off),col="blue" 
          ,density=25,  main= strsplit(input_fp, "_")[[1]][1] )
dev.off()

#axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)




