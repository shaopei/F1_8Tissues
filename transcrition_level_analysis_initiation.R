file_dir="~/Box Sync/BN_IGV/"
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/initiation_CA_nonCA/multipleMaxTSNperTunits/")

organ="BN"
# without RNA-seq
tus=read.table(file = paste(organ,"_allReads_TSS_maxTSNsinProteinCodingTunit_-1-0_mat_patSeq.bed",sep=""), header=FALSE)
# with RNA-seq
tus=read.table(file = paste(organ,"_allReads_TSS_maxTSNsinProteinCodingTunit_-1-0_mat_patSeq_f0.5F0.8gencode.vM25.annotation.gene.bed",sep=""), header=FALSE)

colnames(tus)[13:14] = c("B6_allele", "CAST_allele")
# keep tus with CA in one allele
tus=tus[(tus$B6_allele == "CA" | tus$CAST_allele == "CA"),]
# two group
# One CA, one nonCA
tus$group = "2CA"
tus$nonCA=tus$CAST_allele
tus$nonCA[tus$CAST_allele=="CA"] = tus$B6_allele[tus$CAST_allele=="CA"]
tus$group[(tus$nonCA!= "CA")] = "1CA"
dim(tus)

# keep only -1 is "C"
tus=tus[(tus$nonCA=="CT"|tus$nonCA=="CC"|tus$nonCA=="CG" | tus$nonCA =="CA"),]
dim(tus)
sum(tus$group == "1CA") 
sum(tus$group == "2CA")

# per Tunits, if one of the maxTSN is C{!A}, assign the Tunit to the group 1CA
for (i  in 1:dim(tus)[1]){
  group <- NULL
  if(sum(tus$V10==tus$V10[i]) > 1){
    cat (i)
    cat ("\n")
    if ("1CA" %in% tus$group[tus$V10==tus$V10[i]]){
      tus$group[tus$V10==tus$V10[i]] = "1CA"
    }
  }
}

# per 1CA/2CA group, tunit can only be count once
tus=tus[!duplicated(tus[,c(7:10,which(colnames(tus) == "group"))]),] 

dim(tus)
sum(tus$group == "1CA") 
sum(tus$group == "2CA")


colnames(tus)[8] = "TXSTART"
colnames(tus)[9] = "TXEND"
colnames(tus)[12] = "TXSTRAND"
tus <- tus[(tus$TXEND-tus$TXSTART)>=500,]
dim(tus)
sum(tus$group == "1CA") 
sum(tus$group == "2CA")
# body exclude the first h bp of transcript
h=300
bodies <- tus[,7:12]
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+h
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-h

countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
  pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
  mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))
  
  counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
  if(rpkm==TRUE) {
    counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
  }
  
  return(counts)
}

b6.bw=paste(organ, "_map2ref_1bpbed_map5_B6", sep="")
cast.bw=paste(organ, "_map2ref_1bpbed_map5_CAST", sep="")
filenames <- c(b6.bw, cast.bw)
## Gets counts
counts <- NULL
for(f in filenames) {
  counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE, path=file_dir))
}

tus$B6_counts = counts[,1]
tus$CAST_counts = counts[,2]
tus$log2_CA_NonCA_ratio = log2((tus$B6_counts+1)/(tus$CAST_counts+1))
tus$log2_CA_NonCA_ratio[tus$CAST_allele == "CA"] = log2((tus$CAST_counts+1)/(tus$B6_counts+1))[tus$CAST_allele == "CA"]

library("vioplot")
vioplot(tus$log2_CA_NonCA_ratio ~ tus$group,
        main=paste(organ, " Tunit Proseq", sep=""))
abline(h=0, col="red")
boxplot(tus$log2_CA_NonCA_ratio ~ tus$group,
        #ylim=c(-1,1),
        #xlab="AT long allele",
        ylab="log2(CA+1/NonCA+1)",
        outline=F,las=1,
        main=paste(organ, " Tunit Proseq", sep=""))
abline(h=0, col="red")

dim(tus)
sum(tus$group == "1CA") 
sum(tus$group == "2CA")
wilcox.test(tus$log2_CA_NonCA_ratio[tus$group=="1CA"],
              tus$log2_CA_NonCA_ratio[tus$group=="2CA"])

load(paste("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/RNA-seq/",organ, "_gencode.vM25.annotation.gene-Readcounts.RData", sep=""))
#View(merge.counts)
dim(merge.counts)
colnames(tus)[18]="geneID"
tus_gene = merge(tus , merge.counts, by = "geneID")
dim(tus_gene)
#View(tus_gene)
sum(tus_gene$group == "1CA")
sum(tus_gene$group == "2CA")

tus_gene$log2_CA_NonCA_Exonratio = log2((tus_gene$B6.exon+1)/(tus_gene$CAST.exon+1))
tus_gene$log2_CA_NonCA_Exonratio[tus_gene$CAST_allele == "CA"] = log2((tus_gene$CAST.exon+1)/(tus_gene$B6.exon+1))[tus_gene$CAST_allele == "CA"]

library("vioplot")
vioplot(tus_gene$log2_CA_NonCA_Exonratio~tus_gene$group,
        main=paste(organ, " Gene exon RNA-seq", sep=""))
abline(h=0, col="red")
boxplot(tus_gene$log2_CA_NonCA_Exonratio~tus_gene$group,
        #ylim=c(-1,1),
        #xlab="AT long allele",
        ylab="log2(CA+1/NonCA+1)",
        outline=F,las=1,
        main=paste(organ, " Gene exon RNA-seq", sep=""))
abline(h=0, col="red")

dim(tus_gene)
sum(tus_gene$group == "1CA")
sum(tus_gene$group == "2CA")
wilcox.test(tus_gene$log2_CA_NonCA_Exonratio[tus_gene$group=="1CA"],
            tus_gene$log2_CA_NonCA_Exonratio[tus_gene$group=="2CA"])



boxplot(tus_gene$log2_CA_NonCA_ratio ~ tus_gene$group,
        #ylim=c(-1,1),
        #xlab="AT long allele",
        ylab="log2(CA+1/NonCA+1)",
        outline=F,las=1,
        main=organ)
abline(h=0, col="red")

plot(log2((tus_gene$B6_counts+1)/(tus_gene$CAST_counts+1)), 
  log2((tus_gene$B6.proseq+1)/(tus_gene$CAST.proseq+1)), main=organ)

cor.test(log2((tus_gene$B6_counts+1)/(tus_gene$CAST_counts+1)), 
    log2((tus_gene$B6.proseq+1)/(tus_gene$CAST.proseq+1)))

plot(log2((tus_gene$B6_counts+1)/(tus_gene$CAST_counts+1)), 
     log2((tus_gene$B6.rna+1)/(tus_gene$CAST.rna+1)), main=organ)

plot(log2((tus_gene$B6.proseq+1)/(tus_gene$CAST.proseq+1)), 
     log2((tus_gene$B6.rna+1)/(tus_gene$CAST.rna+1)), 
     main=organ,
     xlim=c(-5,5),
     ylim=c(-7,7))

plot(log(tus_gene$B6.proseq), 
     log(tus_gene$B6.rna), 
     main=organ)

plot(log(tus_gene$CAST.proseq), 
     log(tus_gene$CAST.rna), 
     main=organ)

cor.test(log(tus_gene$CAST.proseq+1), 
     log(tus_gene$CAST.rna+1))


cor.test(log2((tus_gene$B6.proseq+1)/(tus_gene$CAST.proseq+1)), 
         log2((tus_gene$B6.rna+1)/(tus_gene$CAST.rna+1)))

plot(log2((tus_gene$B6.exon+1)/(tus_gene$CAST.exon+1)), 
     log2((tus_gene$B6.rna+1)/(tus_gene$CAST.rna+1)), main=organ)
