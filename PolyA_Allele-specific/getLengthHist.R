#domain_length.pdf
#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/")

args=(commandArgs(TRUE))

folder = args[1]
f1_p = args[2] #"BN_AT_3tunitIntersectNativeHMM_intersectRegion.bed"
f2_p = args[3] #"BN_AT_3tunitIntersectNativeHMM_tunit.bed"
t=args[4]

setwd(folder) #("/Volumes/SPC_SD/IGV/PolyA_Allele-specific")
pdf(paste(t, "length.pdf", sep = "_"))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
f1=read.table(f1_p,header=F)
h<- hist(log10(f1$V3-f1$V2),#col="red" 
#         ,density=25
         , breaks = seq(0,7,1)
#         , freq = T
#         , prob=TRUE
#         ,ylab="Proportion of domains"
#         , xlim=c(0,7)
#         ,xlab="AlleleHMM block length (log10)"
#         ,main= t
#         ,add=F
#         ,las=2
         ,plot =FALSE
         ,right = FALSE
)

h$counts=h$counts/sum(h$counts)

f2=read.table(f2_p,header=F)
h2<- hist(log10(f2$V3-f2$V2),col="blue" 
          , breaks = seq(0,8,1)
          , freq = F
          #,add=F
          #,las=2
          ,plot =FALSE
          ,right = FALSE
)

h2$counts=h2$counts/sum(h2$counts)
plot(h,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,ylim=c(0,max(h$density,h2$density))
     ,xlab="block length"     
     ,las=2
     ,xaxt='n',main= t)
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

plot(h2,col="blue" , add=T)
plot(h,col="red" ,density=25, add=T)

legend("topleft", 
       legend = c( "AT window","Tunit"), 
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25, 10000),
       angle=c(45, 180),
       #angle=45,
       fill=c("red","blue")
       , bty = "n"
)
dev.off()

