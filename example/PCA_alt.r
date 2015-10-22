
#PCA method 1, requires snpStat
require(snpStats)
dat <- read.plink("example/dat/example")
xxmat <- xxt(dat$genotypes, correct.for.missing=F)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors

#PCA method 2, need to have the raw genotype file (plink --bfile example/dat/example --recodeA --out example/result/geno)
dat.geno <- read.table("example/result/geno.raw", header=T)
m <- as.matrix(dat.geno[, -c(1:6)])
for(i in 1:dim(m)[2]){
  mean <- mean(m[,i], na.rm=T) 
  p <- (1 + sum(m[,i], na.rm=T) )/(2+2*sum(!is.na(m[,i])))
  m[,i][is.na(m[,i])] <- mean
  m[,i] <- m[,i] - mean
  m[,i] <- m[,i] / sqrt(p*(1-p))
}
SVD <- svd(m)
pcs <- SVD$u

#read PC from efficientPCA, make sure you have completed the tutorial
pcs_ePCA <- read.table("example/result/pcs.txt", header=F)

#compare PCs
pdf("PCs.pdf", width=10, height=5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(pcs[,1], pcs_ePCA[,1], xlab="PC1 (R)", ylab="PC1 (ePCA)")
plot(pcs[,2], pcs_ePCA[,2], xlab="PC2 (R)", ylab="PC2 (ePCA)")
dev.off()