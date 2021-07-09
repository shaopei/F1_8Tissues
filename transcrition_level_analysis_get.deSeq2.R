##
## deSeq2 ...

require(bigWig)

load("data-counts.RData")

tus <- read.table("../annotations/tuSelecter/final_tus.txt", header=TRUE)
tus <- tus[(tus$TXEND-tus$TXSTART)>500,]
geneID_name <- read.table("../annotations/gencode.vM20_geneID_name_pair.txt", header=F)
#colnames(geneID_name)=c("GENEID", "GENENAME")

library("sqldf")
new_tus <- sqldf ("select TXCHROM,TXSTART,TXEND,GENEID,TXNAME,TXSTRAND,V2,TXTYPE from tus left join geneID_name on tus.GENEID=geneID_name.V1" )
colnames(new_tus)[7]="GENENAME"
tus <- new_tus
bodies <- tus

library("DESeq2")

FDR.lv<- 0.01

## Add specific genes to the plot.
addlab <- function(GENEID, deRes, genes, ...) {  
 idx<-sapply(GENEID, function(GENEID) {
  #ig <- which(genes[,4] == GENEID)
  ig <- which(genes[,7] == GENEID) #use GENENAME instead of GENEID
  io <- ig[which.min(deRes$padj[ig])]
  if(NROW(io)>0) {
        text(deRes$baseMean[io], deRes$log2FoldChange[io], labels= GENEID, cex= 1, pos= 3, ...)
        points(deRes$baseMean[io], deRes$log2FoldChange[io], col="blue", cex=1.5)
        io
  }
 })
 print(idx)
 idx <- unlist(idx); idx <- idx[!is.null(idx)]
 return(data.frame(Gene= genes[idx,7], AveExpr= deRes$baseMean[idx], logFC= deRes$log2FoldChange[idx], adj.P.Val= deRes$padj[idx]))
}

condition <- c(rep("WT", 4), rep("MUT", 4))
replicate <- factor(rep(c(1:4), 2))
Design <- data.frame(colnames(counts), condition, replicate)

## Workaround for deSeq2 issues ...
#setClassUnion("characterORNULL", c("character", "NULL"))
#setClassUnion("DataTableORNULL", c("DataTable", "NULL"))


## Create DESeq2 object.
dds <- DESeqDataSetFromMatrix(countData= counts, colData= Design, design= ~condition+replicate)

dds <- DESeq(dds)
res_wt_mut <- results(dds, contrast=c("condition","WT","MUT"))



 
print(paste("Number of changes: ", sum(res_wt_mut$padj < FDR.lv, na.rm=TRUE))) ## Number of transcripts.
#print(paste("Number of changes: ", sum(res_p_d$padj < FDR.lv, na.rm=TRUE))) ## Number of transcripts.

 pdf(paste("WT-MUT.pdf", sep=""), useDingbats=FALSE)
  plotMA(res_wt_mut, alpha= FDR.lv, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5))
  #addlab(c("ENSMUSG00000028023.16"), res_wt_mut, bodies) ## ENSMUSG00000028023.16 is PITX2
 dev.off()

 pdf(paste("WT-MUT-labeled.pdf", sep=""), useDingbats=FALSE)
  plotMA(res_wt_mut, alpha= FDR.lv, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5))
  addlab(c("Pitx2", "Tbx5"), res_wt_mut, bodies) ## ENSMUSG00000028023.16 is PITX2
  addlab(c("Tbx5"), res_wt_mut, bodies)
  addlab(c("Hcn1"), res_wt_mut, bodies)
  addlab(c("Ryr2"), res_wt_mut, bodies)
  addlab(c("Gja1"), res_wt_mut, bodies)
  #addlab(c("ENSMUSG00000028023.16"), res_wt_mut, bodies) ## ENSMUSG00000028023.16 is PITX2
 dev.off()

 pdf(paste("WT-MUT-fdr0.01-labeled-3.pdf", sep=""), useDingbats=FALSE)
  #plotMA(res_wt_mut[8900:8999,], alpha= FDR.lv, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5))
  #plot(1,1000, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5), xlim = c(0,10000), log="x")
    plotMA(res_wt_mut[c(which(res_wt_mut$padj < FDR.lv),1:10),], alpha= FDR.lv, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5))
    addlab(unique(bodies$GENENAME[res_wt_mut$padj < FDR.lv]), res_wt_mut, bodies)
dev.off()

## Write out genes.
PX <- cbind(bodies, res_wt_mut)
PX <- PX[order(PX$padj),]
write.table(PX, "WT-MUT-changed.genes.tsv", sep="\t", row.names=FALSE, col.names=T, quote=FALSE)
#write.table(PX, "WT-MUT-changed.genes.tsv", sep="\t", row.names=T, col.names=T, quote=FALSE)


