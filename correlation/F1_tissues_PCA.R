setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/correlation")
require(bigWig)
require(cluster)

load("data-rpkms.RData")
load("data-counts.RData")

min_count = 5
indx_counts <- rowSums(counts>=min_count) >= dim(counts)[2]  #every sample has at least 10 reads
indx_trxSize<- (bodies[,3]-bodies[,2])>10000  # to get a robost signal
indx <- indx_counts & indx_trxSize
rpkm <- rpkm[indx,]

# use geneE for cluster


# PCA
pca <- prcomp(rpkm)
library("RColorBrewer")
display.brewer.pal(n = 9, name = "Set1")
scols = brewer.pal(n = 9, name = "Set1")[c(1:5,7:9)]
spch <- c(1,2,3,16,0,4,7,6)

cols<-c()
for (c in scols){
  cols <- c(cols, rep(c,7))
}
# the two problematic lung samples
cols[22]="red"
cols[25]="red"
pch <- c()
for (p in spch){
  pch <- c(pch, rep(p,7))
}

data.frame(colnames(rpkm), cols, pch)
summary(pca)
pc1= "PC1 (52%)"
pc2= "PC2 (18%)"
pc3= "PC3 (8%)"
pc4= "PC4 (5%)"
#pdf("PCA.pdf")
plot(y= pca$rotation[,1], x= pca$rotation[,2], col=cols, pch=pch, xlab=pc2, ylab=pc1)
new_order=c(1,3,6,7,2,8,5,4)
legend("topright",col=scols[new_order],pch=spch[new_order],
       legend=stage[new_order], bty="n")
pairs(pca$rotation[,1:5], col=cols, pch=pch)

plot(y= pca$rotation[,1], x= pca$rotation[,4], col=cols, pch=pch, xlab=pc4, ylab=pc1)
legend("bottomright",col=scols[new_order],pch=spch[new_order],
       legend=stage[new_order], bty="n")

plot(y= pca$rotation[,2], x= pca$rotation[,3], col=cols, pch=pch, xlab=pc3, ylab=pc2)
legend("topright",col=scols[new_order],pch=spch[new_order],
       legend=stage[new_order], bty="n")


# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
mydata <- t(rpkm)
d <- dist(mydata) # euclidean distances between the rows
fit <- cmdscale(d, eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS", type="p", col=cols, pch=pch)
legend("topright",col=scols[new_order],pch=spch[new_order],
       legend=stage[new_order], bty="n")


# PCA exclude SK
pca <- prcomp(rpkm[c(1:35,43:56),])

summary(pca)
#pdf("PCA.pdf")
cols=cols[c(1:35,43:56)]
pch=pch[c(1:35,43:56)]

plot(y= pca$rotation[,1], x= pca$rotation[,2], col=cols, pch=pch, xlab="PC2", ylab="PC1")
legend("topleft",col=scols[new_order],pch=spch[new_order],
       legend=stage[new_order], bty="n")

pairs(pca$rotation[,1:5], col=cols, pch=pch)

plot(y= pca$rotation[,1], x= pca$rotation[,4], col=cols, pch=pch, xlab="PC4", ylab="PC1")
legend("topleft",col=scols[new_order],pch=spch[new_order],
       legend=stage[new_order], bty="n")



