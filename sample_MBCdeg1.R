#########################################
### This file can be used to execute MBCdeg1.
#########################################
set.seed(1)
###  Parameters  ###
in_f <- "sample.txt"             #input filename
out_f <- "sample_MBCdeg1.txt"    #output filename
n1 <- 5                          #number of replicates for group 1 (or A)
n2 <- 6                          #number of replicates for group 2 (or B)
K <- 3                           #preselected number of clusters

###  Load packages  ###
library(MBCluster.Seq)

###  Read data  ###
data <- read.table(in_f, header=TRUE, row.names=1, sep="\t")
dim(data)
data.cl <- c(rep(1, n1), rep(2, n2))

###  MBCdeg1 (main)  ###
hoge <- RNASeq.Data(data, Normalizer=NULL,
           Treatment=data.cl, GeneID=rownames(data))
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")

###  Identify the k value that corresponds to the non-DEG cluster  ###
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)

###  Output (file)  ###
pattern <- NULL
pattern[1:K] <- "DEG2"
pattern[cls$centers[,1] > 0] <- "DEG1"
pattern[k] <- "non-DEG"
pat <- cls$cluster
for(i in 1:K){
  pat[pat == i] <- pattern[i]
}
PP <- cls$probability
PP_nonDEG <- cls$probability[,k]
ranking <- rank(cls$probability[,k])
clust_num <- cls$cluster
tmp <- cbind(rownames(data), data, PP, PP_nonDEG, ranking, pat, clust_num)
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output (R console)  ###
table(pat)
table(clust_num)
cls$centers
