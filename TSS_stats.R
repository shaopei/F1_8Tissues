setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSS/")
t="BN"
tss_dir="./"
name_body="_allReads_TSS"
TSS <- read.table(paste(tss_dir,t,name_body,".bed", sep =""), header = F)
names(TSS)[4]="TSNCount"
names(TSS)[5]="TSSReadsCount"


plot(log10(TSS$TSSReadsCount),TSS$TSNCount,  
     xlim=c(0,4), 
     ylim=c(0,80), 
     pch=19, col=rgb(0,0,0,alpha = 0.125/2))


library("ggplot2")

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
myColor_scale_fill <- scale_fill_gradientn(colours = myColor, trans='log10')

p <- ggplot(TSS, aes(x=log10(TSSReadsCount), y=TSNCount)) +
  geom_bin2d(bins = 150) +  myColor_scale_fill + theme_bw() +
  xlim(1, 4) + ylim(0,180) + labs(x="log10 (TSS Read Counts)")  + labs(y="TSN counts")
p
p +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=20), #,face="bold"))
  legend.title = element_text(size = 20),
  legend.text = element_text(size = 14)
)

ggsave('~/Desktop/temp.pdf', p, h=7, w=7)