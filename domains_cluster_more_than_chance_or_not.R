setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/domains_cluster_more_than_chance_or_not_tunit_protein/")
Head="LV"
df_0.1=read.table(paste(Head, "_TunitProteinSrainEffect_binomtest_fdr0.1_adjacentTunit.bed", sep=""))
df_0.9=read.table(paste(Head,"_TunitProteinSrainEffect_binomtest_fdr0.9_adjacentTunit.bed", sep=""))

colnames(df_0.1)[22]="pair_fdr"  # fdr of the tunit that pair with the fdr0.1 tunit
colnames(df_0.1)[23]="pair_distance" # distance of the tunit that pair with the fdr0.1 tunit
df_0.1$TunitLength = df_0.1$V3 - df_0.1$V2
df_0.1$expLevel = df_0.1$V18 + df_0.1$V19 + df_0.1$V20

colnames(df_0.9)[22]="pair_fdr"
colnames(df_0.9)[23]="pair_distance"
df_0.9$TunitLength = df_0.9$V3 - df_0.9$V2
df_0.9$expLevel = df_0.9$V18 + df_0.9$V19 + df_0.9$V20


#Do TUs that are adjacent to a significant change tend to share the same change?
fisher.test(matrix(c(dim(df_0.9)[1],sum(df_0.9$pair_fdr <= 0.1), 
                    dim(df_0.1)[1],sum(df_0.1$pair_fdr <= 0.1)
                     ) , 2,2))


# Do TUs change in the same direction of change? Or do they often have the opposite direction?
df = df_0.1[df_0.1$pair_fdr <= 0.1,]
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V4), ",")))
df$Tunit1_winP = temp[,1]
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V15), ",")))
df$Tunit2_winP = temp[,1]

# same direction
levels(df$Tunit2_winP) = levels(df$Tunit1_winP)
dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1] / dim(df)[1]
# opposite direction
dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1] / dim(df)[1]
pie(c(dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1], dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1] ),
    label=c("same direction", "opposite direction"))
# same direction
dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1]
# opposite direction
dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1]


# Distance distribution
#pdf("Distance_distribution_Between_Tunit_pairs.pdf", width=5, height = 5, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)

h1<- hist(log10(df_0.1$pair_distance[df_0.1$pair_fdr <= 0.1]), plot = F)
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="distance between pair Tunits"     
     ,las=2
     ,xaxt='n',main= ""
     )
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(df_0.9$pair_distance[df_0.9$pair_fdr <= 0.1])
          ,plot =FALSE
)
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "FDR <=0.1","FDR > 0.9"), 
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
ks.test(df_0.9$pair_distance[df_0.9$pair_fdr <= 0.1], 
        df_0.1$pair_distance[df_0.1$pair_fdr <= 0.1] )


# Tunit length
h1<- hist(log10(df_0.1$TunitLength)[!duplicated(df_0.9[,1:6])])
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="Tunit legnth"     
     ,las=2
     ,xaxt='n',main= ""
)
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(df_0.9$TunitLength)[!duplicated(df_0.1[,1:6])]
          ,plot =FALSE
)
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "FDR <=0.1","FDR > 0.9"), 
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

ks.test(df_0.9$TunitLength[!duplicated(df_0.9[,1:6])], 
        df_0.1$TunitLength[!duplicated(df_0.1[,1:6])])

# Tunit expression level
h1<- hist(log10(df_0.1$expLevel)[!duplicated(df_0.9[,1:6])])
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="Tunit Expression Level"     
     ,las=2
     ,xaxt='n',main= ""
)
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(df_0.9$expLevel)[!duplicated(df_0.1[,1:6])]
          ,plot =FALSE
)
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "FDR <=0.1","FDR > 0.9"), 
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

ks.test(df_0.9$expLevel[!duplicated(df_0.9[,1:6])], 
        df_0.1$expLevel[!duplicated(df_0.1[,1:6])])


#####
# two analysis, one use all tunits that are nearest degrardless of strand,
# the other only use nearest from the opposite strand
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/domains_cluster_more_than_chance_or_not_tunit_protein/")

# with AT window VS withOUT AT window
Head="LV"
df_withAT=read.table(paste(Head, "_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunit.bed", sep=""))
df_withoutAT=read.table(paste(Head, "_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_adjacentTunit.bed", sep=""))

#use nearest from the opposite strand
df_withAT=read.table(paste(Head, "_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunitOppositeStrand.bed", sep=""))
df_withoutAT=read.table(paste(Head, "_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_adjacentTunitOppositeStrand.bed", sep=""))
#df_withoutAT = new_df
dim(df_withAT)
df_withAT = unique(df_withAT[,1:23])
dim(df_withAT)
colnames(df_withAT)[22]="pair_fdr"  # fdr of the tunit that pair with the fdr0.1 tunit
colnames(df_withAT)[23]="pair_distance" # distance of the tunit that pair with the fdr0.1 tunit
df_withAT$TunitLength = df_withAT$V3 - df_withAT$V2
df_withAT$expLevel = df_withAT$V18 + df_withAT$V19 + df_withAT$V20

colnames(df_withoutAT)[22]="pair_fdr"
colnames(df_withoutAT)[23]="pair_distance"
df_withoutAT$TunitLength = df_withoutAT$V3 - df_withoutAT$V2
df_withoutAT$expLevel = df_withoutAT$V18 + df_withoutAT$V19 + df_withoutAT$V20

# go to current+179=368 line # subsample withoutAT to match the distrubution of expLevel in withAT 
#df_withoutAT = new_df 

resample=TRUE
set.seed(100)

if(resample){
    # distribution of expLevel withAT
    # set the breaks to the same length
    h1<- hist(log10(df_withAT$expLevel)[!duplicated(df_withAT[,1:6])], breaks = seq(0,7,0.5))
    h1$frac=h1$counts/sum(h1$counts)
    # distribution of expLevel withoutAT
    h2<- hist(log10(df_withoutAT$expLevel)[!duplicated(df_withoutAT[,1:6])]
              , breaks = seq(0,7,0.5))
    
    
    # use withAT fraction to calculate the counts of tunits to be sampled from each expLevel breaks
    # determine *sum(h2$counts)/10 by observing numbers in excel
    f = round(h1$frac*sum(h2$counts)/10) 
    # the breaks are not the same length
    # h2$counts = h2$counts[(1+length(h2$counts) - length(f)):length(h2$counts)] 
    
    #f[f > h2$counts] = h2$counts[f > h2$counts]
    
    df_withoutAT$log10expLevel = log10(df_withoutAT$expLevel)
    #u=1.5
    #l=1
    new_df <- NULL
    
    for (i in 1:(length(h1$breaks)-1)){
        l = h1$breaks[i]
        u = h1$breaks[i+1]
        temp=df_withoutAT[df_withoutAT$log10expLevel>l & df_withoutAT$log10expLevel<=u,]
        cat ("draw",f[i],"from" ,h2$counts[i], "\n")
        if (f[i] <= h2$counts[i]){
            t1 = temp[sample(nrow(temp),f[i], replace = FALSE), ]
        }else{
            t1 = temp[sample(nrow(temp),f[i], replace = TRUE), ]
        }
        new_df=rbind.data.frame(new_df, t1)
    }
    #View(new_df)
    
    h1<- hist(log10(df_withAT$expLevel)[!duplicated(df_withAT[,1:6])])
    h1$counts=h1$counts/sum(h1$counts)
    plot(h1,col="red" 
         ,density=25     
         ,ylab="Proportion"
         , xlim=c(0,7)
         ,xlab="log10(Tunit Expression Level)"     
         ,las=2
         #,xaxt='n',main= ""
    )
    #axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)
    
    h2<- hist(log10(new_df$expLevel)[!duplicated(new_df[,1:6])]
              ,plot =FALSE
    )
    h2$counts=h2$counts/sum(h2$counts)
    plot(h2,col="blue" , add=T)
    plot(h1,col="red" ,density=25, add=T)
    
    legend("topright", 
           legend = c( "With AT window","Without AT window"), 
           #pch=c(15,15),
           #cex=2, 
           lty=c(0,0),
           #bty="n",
           lwd=1.5, 
           density=c(25, 10000),
           angle=c(45, 180),
           #angle=45,
           fill=c("red","blue")
           , bty = "n"
    )
    
    ks.test(new_df$expLevel[!duplicated(new_df[,1:6])], 
            df_withAT$expLevel[!duplicated(df_withAT[,1:6])])
    df_withoutAT = new_df 
}

#Given a gene with AT window,  is the adjacent gene more likely to be biased? 
fisher.test(matrix(c(dim(unique(df_withoutAT[,1:6]))[1],dim(unique(df_withoutAT[df_withoutAT$pair_fdr <= 0.1, 1:6]))[1], 
                     dim(unique(df_withAT[,1:6]))[1],dim(unique(df_withAT[df_withAT$pair_fdr <= 0.1, 1:6]))[1])
                   , 2,2))


# Do TUs change in the same direction of change? Or do they often have the opposite direction?
df = df_withAT[df_withAT$pair_fdr <= 0.1,]
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V4), ",")))
df$Tunit1_winP = temp[,1]
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V15), ",")))
df$Tunit2_winP = temp[,1]

# same direction
levels(df$Tunit2_winP) = levels(df$Tunit1_winP)
dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1] / dim(df)[1]
# opposite direction
dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1] / dim(df)[1]
pie(c(dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1], dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1] ),
    label=c("same direction", "opposite direction"))
# same direction
dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1]
# opposite direction
dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1]


# Distance distribution
#pdf("Distance_distribution_Between_Tunit_pairs_withVSwithoutATwindow.pdf", width=5, height = 5, useDingbats=FALSE)
#par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
#par(mgp=c(3,1,0))
#par(cex.lab=2.2, cex.axis=2.2)

h1<- hist(log10(df_withAT$pair_distance[df_withAT$pair_fdr <= 0.1]),breaks=seq(0,7,0.5))
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="distance between pair Tunits"     
     ,las=2
     ,xaxt='n',main= ""
     #,breaks=seq(1,7,0.5)
     )

axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(df_withoutAT$pair_distance[df_withoutAT$pair_fdr <= 0.1])
          ,breaks=seq(0,7,0.5)
          #,plot =FALSE
)
h2$counts=h2$counts/sum(h2$counts)
plot(h1,col="red" ,density=25)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "With AT window","Without AT window"), 
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

#dev.off()
ks.test(df_withoutAT$pair_distance[df_withoutAT$pair_fdr <= 0.1], 
        df_withAT$pair_distance[df_withAT$pair_fdr <= 0.1] )


# Tunit length
h1<- hist(log10(df_withAT$TunitLength)[!duplicated(df_withAT[,1:6])], 
          breaks=seq(2,7,0.5))
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="Tunit legnth"     
     ,las=2
     ,xaxt='n',main= ""
)
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(df_withoutAT$TunitLength)[!duplicated(df_withoutAT[,1:6])]
          ,plot =FALSE,breaks=seq(2,7,0.5)
)
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "With AT window","Without AT window"), 
       #pch=c(15,15),
#       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25, 10000),
       angle=c(45, 180),
       #angle=45,
       fill=c("red","blue")
       , bty = "n"
)

ks.test(df_withoutAT$TunitLength[!duplicated(df_withoutAT[,1:6])], 
        df_withAT$TunitLength[!duplicated(df_withAT[,1:6])])

# Tunit expression level
h1<- hist(log10(df_withAT$expLevel)[!duplicated(df_withAT[,1:6])])
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="Tunit Expression Level"     
     ,las=2
     ,xaxt='n',main= ""
)
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(df_withoutAT$expLevel)[!duplicated(df_withoutAT[,1:6])]
          ,plot =FALSE
)
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "With AT window","Without AT window"), 
       #pch=c(15,15),
    #   cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25, 10000),
       angle=c(45, 180),
       #angle=45,
       fill=c("red","blue")
       , bty = "n"
)

ks.test(df_withoutAT$expLevel[!duplicated(df_withoutAT[,1:6])], 
        df_withAT$expLevel[!duplicated(df_withAT[,1:6])])

# examine counts
h1<- hist(log10(df_withAT$expLevel)[!duplicated(df_withAT[,1:6])])
#h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="Tunit Expression Level"     
     ,las=2
     ,xaxt='n',main= ""
)
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(df_withoutAT$expLevel)[!duplicated(df_withoutAT[,1:6])]
          ,plot =FALSE
)
#h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" )
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "With AT window","Without AT window"), 
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

# subsample withoutAT to match the distrubution of expLevel in withAT
# sample without replacement to avid duplicates 

# distribution of expLevel withAT
# set the breaks to the same length
h1<- hist(log10(df_withAT$expLevel)[!duplicated(df_withAT[,1:6])], breaks = seq(0,7,0.5))
h1$frac=h1$counts/sum(h1$counts)
# distribution of expLevel withoutAT
h2<- hist(log10(df_withoutAT$expLevel)[!duplicated(df_withoutAT[,1:6])]
          , breaks = seq(0,7,0.5))


# use withAT fraction to calculate the counts of tunits to be sampled from each expLevel breaks
# determine *sum(h2$counts)/10 by observing numbers in excel
f = round(h1$frac*sum(h2$counts)/10) 
# the breaks are not the same length
# h2$counts = h2$counts[(1+length(h2$counts) - length(f)):length(h2$counts)] 

#f[f > h2$counts] = h2$counts[f > h2$counts]

df_withoutAT$log10expLevel = log10(df_withoutAT$expLevel)
#u=1.5
#l=1
new_df <- NULL

for (i in 1:(length(h1$breaks)-1)){
    l = h1$breaks[i]
    u = h1$breaks[i+1]
    temp=df_withoutAT[df_withoutAT$log10expLevel>l & df_withoutAT$log10expLevel<=u,]
    cat ("draw",f[i],"from" ,h2$counts[i], "\n")
    if (f[i] <= h2$counts[i]){
        t1 = temp[sample(nrow(temp),f[i], replace = FALSE), ]
    }else{
        t1 = temp[sample(nrow(temp),f[i], replace = TRUE), ]
    }
    new_df=rbind.data.frame(new_df, t1)
}
View(new_df)

h1<- hist(log10(df_withAT$expLevel)[!duplicated(df_withAT[,1:6])])
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,xlab="log10(Tunit Expression Level)"     
     ,las=2
     #,xaxt='n',main= ""
)
#axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

h2<- hist(log10(new_df$expLevel)[!duplicated(new_df[,1:6])]
          ,plot =FALSE
)
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "With AT window","Without AT window"), 
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25, 10000),
       angle=c(45, 180),
       #angle=45,
       fill=c("red","blue")
       , bty = "n"
)

ks.test(new_df$expLevel[!duplicated(new_df[,1:6])], 
        df_withAT$expLevel[!duplicated(df_withAT[,1:6])])

### compare genes (use protein coding Tunits here) overlap with AT window, and genes overlap with gene(another tunit) without AT window
# Do TUs that are TSS within AT window of another TU more likely to be biased? 
# with AT window VS withOUT AT window
Head="BN"
df_withAT=read.table(paste(Head, "_TunitProteinSrainEffect_binomtest_fdrAll_withATwindow_adjacentTunitOppositeStrand.bed", sep=""))
df_withoutAT=read.table(paste(Head, "_TunitProteinSrainEffect_binomtest_fdrAll_withoutATwindow_adjacentTunitOppositeStrand.bed", sep=""))

colnames(df_withAT)[22]="pair_fdr"  # fdr of the tunit that pair with the fdr0.1 tunit
colnames(df_withAT)[23]="pair_distance" # distance of the tunit that pair with the fdr0.1 tunit
df_withAT$TunitLength = df_withAT$V3 - df_withAT$V2
df_withAT$expLevel = df_withAT$V7 + df_withAT$V8 + df_withAT$V9

colnames(df_withoutAT)[22]="pair_fdr"
colnames(df_withoutAT)[23]="pair_distance"
df_withoutAT$TunitLength = df_withoutAT$V3 - df_withoutAT$V2
df_withoutAT$expLevel = df_withoutAT$V7 + df_withoutAT$V8 + df_withoutAT$V9

# plus strand
df_withoutAT$TSSin2ndPart[df_withoutAT$V6 =="+"] = df_withoutAT$V14[df_withoutAT$V6 =="+"] <= df_withoutAT$V3[df_withoutAT$V6 =="+"] &  df_withoutAT$V14[df_withoutAT$V6 =="+"] >= ((df_withoutAT$V3+df_withoutAT$V2)/2)[df_withoutAT$V6 =="+"]

# minus strand 
df_withoutAT$TSSin2ndPart[df_withoutAT$V6 =="-"] = df_withoutAT$V13[df_withoutAT$V6 =="-"] >= df_withoutAT$V2[df_withoutAT$V6 =="-"] &  df_withoutAT$V13[df_withoutAT$V6 =="-"] <= ((df_withoutAT$V3+df_withoutAT$V2)/2)[df_withoutAT$V6 =="-"]
dim(df_withoutAT)
df_withoutAT = df_withoutAT[df_withoutAT$TSSin2ndPart,]
dim(df_withoutAT)

#Do TUs that are within AT window more likely to be biased?
fisher.test(matrix(c(dim(unique(df_withoutAT[,1:6]))[1],dim(unique(df_withoutAT[df_withoutAT$pair_fdr <= 0.1, 1:6]))[1], 
                     dim(unique(df_withAT[,1:6]))[1],dim(unique(df_withAT[df_withAT$pair_fdr <= 0.1, 1:6]))[1])
                   , 2,2))

# Do TUs change in the same direction of change? Or do they often have the opposite direction?
df = df_withAT[df_withAT$pair_fdr <= 0.1,]
df = df_withoutAT[df_withoutAT$pair_fdr <= 0.1,]
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V4), ",")))
df$Tunit1_winP = temp[,1]
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V15), ",")))
df$Tunit2_winP = temp[,1]

# same direction
levels(df$Tunit2_winP) = levels(df$Tunit1_winP)
dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1] / dim(unique(df[,1:6]))[1]
# opposite direction
dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1] / dim(unique(df[,1:6]))[1]
pie(c(dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1], dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1] ),
    label=c("same direction", "opposite direction"))
# same direction
dim(unique(df[df$Tunit1_winP == df$Tunit2_winP, 1:6]))[1]
# opposite direction
dim(unique(df[df$Tunit1_winP != df$Tunit2_winP, 1:6]))[1]

fisher.test(matrix(c(72+68,68, 
                     89+81,81)
                   , 2,2))





### gene annotation counts OR tunit counts in Strain effect domain with VS without AT window
Head="BN"
df_with = read.table(paste(Head, "_T8_2Strand_p0.05_effect_strain.bed_cluster_withATwindow_geneCount.txt", sep=""))
df_without = read.table(paste(Head, "_T8_2Strand_p0.05_effect_strain.bed_cluster_withoutATwindow_geneCount.txt", sep=""))

h1<- hist(df_with$V1, breaks = seq(0,200,5))
h1$counts=h1$counts/sum(h1$counts)
plot(h1,col="red" 
     ,density=25     
     ,ylab="Proportion"
     #, xlim=c(0,7)
     #,xlab="Tunit Expression Level"     
     ,las=1
     #,xaxt='n'
     ,main= Head
)

h2<- hist(df_without$V1
          ,plot =FALSE
          , breaks = seq(0,200,5))
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h1,col="red" ,density=25, add=T)

legend("topright", 
       legend = c( "With AT window","Without AT window"), 
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25, 10000),
       angle=c(45, 180),
       #angle=45,
       fill=c("red","blue")
       , bty = "n"
)

ks.test(df_without$V1, df_with$V1)
        

        
        
        