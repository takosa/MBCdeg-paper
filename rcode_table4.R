#########################################
### This file can be used to obtain the 
### similar results in Table 4 and Additional file 8.
#########################################
set.seed(1)
###  Parameters  ###
param_ID <- "SRP001540"          #ID for Pickrell's data
n1 <- 40                         #number of replicates for group 1
n2 <- 29                         #number of replicates for group 2
FDR <- 0.1                       #FDR threshold
out1 <- "Table4.txt"             #output filename 1
out2 <- "Additional8_Sheet1.txt" #output filename 2
out3 <- "Additional8_Sheet2.txt" #output filename 3

###  Load packages  ###
library(edgeR)
library(DESeq2)
library(TCC)
library(MBCluster.Seq)
library(recount)

###  Download and read the data  ###
download_study(param_ID, type="rse-gene", download=T)
load(file.path(param_ID, 'rse_gene.Rdata'))
rse <- scale_counts(rse_gene)
x <- assays(rse)$counts
x <- as.data.frame(x)

###  Collapsing the data for technical replicates  ###
data <- cbind(x$SRR031822, x$SRR031953 + x$SRR031873,#Female1-2
              x$SRR031952 + x$SRR031871, x$SRR031868,#Female3-4
              x$SRR031819, x$SRR031897 + x$SRR031857,#Female5-6
              x$SRR031823, x$SRR031959, x$SRR031955,#Female7-9
              x$SRR031954, x$SRR031956, x$SRR031838,#Female10-12
              x$SRR031918, x$SRR031817, x$SRR031949 + x$SRR031852,#Female13-15
              x$SRR031841, x$SRR031865, x$SRR031896,#Female16-18
              x$SRR031853, x$SRR031820, x$SRR031874,#Female19-21
              x$SRR031895, x$SRR031870, x$SRR031839,#Female22-24
              x$SRR031958, x$SRR031867, x$SRR031848,#Female25-27
              x$SRR031847, x$SRR031818, x$SRR031919,#Female28-30
              x$SRR031866, x$SRR031849, x$SRR031877,#Female31-33
              x$SRR031814, x$SRR031914, x$SRR031812,#Female34-36
              x$SRR031842, x$SRR031843, x$SRR031860, x$SRR031837,#Female37-40
              x$SRR031917, x$SRR031821 + x$SRR031898,#Male1-2
              x$SRR031950 + x$SRR031850, x$SRR031876 + x$SRR031862,#Male3-4
              x$SRR031875, x$SRR031915, x$SRR031878 + x$SRR031863,#Male5-7
              x$SRR031869, x$SRR031864, x$SRR031845,#Male8-10
              x$SRR031951 + x$SRR031851, x$SRR031846,#Male11-12
              x$SRR031916, x$SRR031844, x$SRR031813,#Male13-15
              x$SRR031894, x$SRR031854, x$SRR031858,#Male16-18
              x$SRR031859, x$SRR031872, x$SRR031816,#Male19-21
              x$SRR031815, x$SRR031920 + x$SRR031899,#Male22-23
              x$SRR031957 + x$SRR031855, x$SRR031840,#Male24-25
              x$SRR031948, x$SRR031893, x$SRR031811, x$SRR031861)#Male26-29
colnames(data) <- c(paste("Female", 1:40, sep=""), paste("Male", 1:29, sep=""))
rownames(data) <- rownames(x)
dim(data)

###  Filtering low count genes  ###
obj <- as.logical(rowSums(data) > 0)
data <- data[obj,]
dim(data)
counts <- data

###  Preparation for class labels  ###
data.cl <- c(rep(1, n1), rep(2, n2))
group <- data.cl

###  edgeR  ###
y <- edgeR::DGEList(counts = counts, group = group)
y <- edgeR::calcNormFactors(y)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
fit <- edgeR::glmQLFit(y, design)
qlf <- edgeR::glmQLFTest(fit, coef = 2)
res <- edgeR::topTags(qlf, n = nrow(counts), sort.by = "none")
p.value <- res$table$PValue
q.value <- p.adjust(p.value, method="BH")
obj <- (q.value < FDR)
pat <- NULL
pat[1:nrow(data)] <- "DEG2"
pat[res$table$logFC < 0] <- "DEG1"
pat[!obj] <- "non-DEG"
#table(pat)
pat_edger <- pat
logratio_edger <- res$table$logFC

###  DESeq2  ###
colData <- data.frame(condition=as.factor(data.cl))
d <- DESeqDataSetFromMatrix(countData=data,
             colData=colData, design=~condition)
d <- DESeq(d)
tmp <- results(d)
p.value <- tmp$pvalue
p.value[is.na(p.value)] <- 1
q.value <- tmp$padj
q.value[is.na(q.value)] <- 1
obj <- (q.value < FDR)
pat <- NULL
pat[1:nrow(data)] <- "DEG2"
pat[res$table$logFC < 0] <- "DEG1"
pat[!obj] <- "non-DEG"
#table(pat)
pat_deseq2 <- pat

###  TCC  ###
tcc <- new("TCC", data, data.cl)
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=0.05)
result <- getResult(tcc, sort=FALSE)
q.value <- tcc$stat$q.value
obj <- (q.value < FDR)
pat <- NULL
pat[1:nrow(data)] <- "DEG2"
pat[res$table$logFC < 0] <- "DEG1"
pat[!obj] <- "non-DEG"
#table(pat)
pat_tcc <- pat

#Obtaining size factors from TCC's normalization factors
norm.factors <- tcc$norm.factors
ef.libsizes <- colSums(data)*norm.factors
size.factors <- ef.libsizes/mean(ef.libsizes)

###  MBCdeg1 (preparation)  ###
hoge <- RNASeq.Data(data, Normalizer=NULL,
           Treatment=data.cl, GeneID=rownames(data))

###  MBCdeg1 (K=3)  ###
K <- 3
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)
pattern <- NULL
pattern[1:K] <- "DEG2"
pattern[cls$centers[,1] > 0] <- "DEG1"
pattern[k] <- "non-DEG"
pat <- cls$cluster
for(i in 1:K){
  pat[pat == i] <- pattern[i]
}
pat_MBCdeg1_3 <- pat
tmp1 <- paste("cluster", 1:K, sep="")
tmp2 <- paste("mu", 1:K, sep="")
result <- cbind(tmp1, table(cls$cluster), tmp2, cls$centers, L2norm, pattern)
colnames(result) <- c("(a) MBCdeg1(K=3)", "#Genes", "Cluster centers",
                      "Group 1", "Group 2", "L2 Norm", "Pattern")
result_MBCdeg1_3 <- result

###  MBCdeg1 (K=4)  ###
K <- 4
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)
pattern <- NULL
pattern[1:K] <- "DEG2"
pattern[cls$centers[,1] > 0] <- "DEG1"
pattern[k] <- "non-DEG"
pat <- cls$cluster
for(i in 1:K){
  pat[pat == i] <- pattern[i]
}
pat_MBCdeg1_4 <- pat
tmp1 <- paste("cluster", 1:K, sep="")
tmp2 <- paste("mu", 1:K, sep="")
result <- cbind(tmp1, table(cls$cluster), tmp2, cls$centers, L2norm, pattern)
colnames(result) <- c("(b) MBCdeg1(K=4)", "#Genes", "Cluster centers",
                      "Group 1", "Group 2", "L2 Norm", "Pattern")
result_MBCdeg1_4 <- result

###  MBCdeg1 (K=5)  ###
K <- 5
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)
pattern <- NULL
pattern[1:K] <- "DEG2"
pattern[cls$centers[,1] > 0] <- "DEG1"
pattern[k] <- "non-DEG"
pat <- cls$cluster
for(i in 1:K){
  pat[pat == i] <- pattern[i]
}
pat_MBCdeg1_5 <- pat
tmp1 <- paste("cluster", 1:K, sep="")
tmp2 <- paste("mu", 1:K, sep="")
result <- cbind(tmp1, table(cls$cluster), tmp2, cls$centers, L2norm, pattern)
colnames(result) <- c("(c) MBCdeg1(K=5)", "#Genes", "Cluster centers",
                      "Group 1", "Group 2", "L2 Norm", "Pattern")
result_MBCdeg1_5 <- result

###  MBCdeg2 (preparation)  ###
hoge <- RNASeq.Data(data, Normalizer=log2(size.factors),
           Treatment=data.cl, GeneID=rownames(data))

###  MBCdeg2 (K=3)  ###
K <- 3
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)
pattern <- NULL
pattern[1:K] <- "DEG2"
pattern[cls$centers[,1] > 0] <- "DEG1"
pattern[k] <- "non-DEG"
pat <- cls$cluster
for(i in 1:K){
  pat[pat == i] <- pattern[i]
}
pat_MBCdeg2_3 <- pat
tmp1 <- paste("cluster", 1:K, sep="")
tmp2 <- paste("mu", 1:K, sep="")
result <- cbind(tmp1, table(cls$cluster), tmp2, cls$centers, L2norm, pattern)
colnames(result) <- c("(d) MBCdeg2(K=3)", "#Genes", "Cluster centers",
                      "Group 1", "Group 2", "L2 Norm", "Pattern")
result_MBCdeg2_3 <- result

###  MBCdeg2 (K=4)  ###
K <- 4
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)
pattern <- NULL
pattern[1:K] <- "DEG2"
pattern[cls$centers[,1] > 0] <- "DEG1"
pattern[k] <- "non-DEG"
pat <- cls$cluster
for(i in 1:K){
  pat[pat == i] <- pattern[i]
}
pat_MBCdeg2_4 <- pat
tmp1 <- paste("cluster", 1:K, sep="")
tmp2 <- paste("mu", 1:K, sep="")
result <- cbind(tmp1, table(cls$cluster), tmp2, cls$centers, L2norm, pattern)
colnames(result) <- c("(e) MBCdeg2(K=4)", "#Genes", "Cluster centers",
                      "Group 1", "Group 2", "L2 Norm", "Pattern")
result_MBCdeg2_4 <- result

###  MBCdeg2 (K=5)  ###
K <- 5
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)
pattern <- NULL
pattern[1:K] <- "DEG2"
pattern[cls$centers[,1] > 0] <- "DEG1"
pattern[k] <- "non-DEG"
pat <- cls$cluster
for(i in 1:K){
  pat[pat == i] <- pattern[i]
}
pat_MBCdeg2_5 <- pat
tmp1 <- paste("cluster", 1:K, sep="")
tmp2 <- paste("mu", 1:K, sep="")
result <- cbind(tmp1, table(cls$cluster), tmp2, cls$centers, L2norm, pattern)
colnames(result) <- c("(f) MBCdeg2(K=5)", "#Genes", "Cluster centers",
                      "Group 1", "Group 2", "L2 Norm", "Pattern")
result_MBCdeg2_5 <- result

###  concatenate the results  ###
hoge <- cbind(pat_edger, pat_deseq2, pat_tcc,
              pat_MBCdeg1_3, pat_MBCdeg1_4, pat_MBCdeg1_5,
              pat_MBCdeg2_3, pat_MBCdeg2_4, pat_MBCdeg2_5)

#################
###  Table 4  ###
#################
#numbers of genes in individual patterns
result <- NULL
for(i in 1:ncol(hoge)){
  result <- rbind(result, table(hoge[,i]))
}
result1 <- result

#mean logratios for individual patterns
result <- NULL
for(i in 1:ncol(hoge)){
  tmp <- c(mean(logratio_edger[hoge[,i] == "DEG1"]),
           mean(logratio_edger[hoge[,i] == "DEG2"]),
           mean(logratio_edger[hoge[,i] == "non-DEG"]))
  result <- rbind(result, tmp)
}
result2 <- result

result <- matrix(0, nrow=nrow(result1), ncol=ncol(result1))
for(i in 1:nrow(result1)){
  for(j in 1:ncol(result1)){
    result[i,j] <- paste(result1[i,j], "(", sprintf("%5.3f", result2[i,j]), ")", sep="")
  }
}
colnames(result) <- c("DEG1", "DEG2", "non-DEG")
Method <- c("edgeR", "DESeq2", "TCC",
                      "MBCdeg1(K=3)", "MBCdeg1(K=4)", "MBCdeg1(K=5)",
                      "MBCdeg2(K=3)", "MBCdeg2(K=4)", "MBCdeg2(K=5)")
tmp <- cbind(Method, result)
write.table(tmp, out1, sep="\t", append=F, quote=F, row.names=F)

####################################
###  Sheet1 in Additional8.xlsx  ###
####################################
result <- matrix(10, nrow=ncol(hoge), ncol=ncol(hoge))
for(i in 1:ncol(hoge)){
  for(j in 1:ncol(hoge)){
    result[i,j] <- sum(hoge[,i] == hoge[,j])/nrow(hoge)
  }
}
colnames(result) <- Method
tmp <- cbind(Method, result)
write.table(tmp, out2, sep="\t", append=F, quote=F, row.names=F)

####################################
###  Sheet2 in Additional8.xlsx  ###
####################################
tmp <- rbind(colnames(result_MBCdeg1_3), result_MBCdeg1_3,
                colnames(result_MBCdeg1_4), result_MBCdeg1_4,
                colnames(result_MBCdeg1_5), result_MBCdeg1_5,
                colnames(result_MBCdeg2_3), result_MBCdeg2_3,
                colnames(result_MBCdeg2_4), result_MBCdeg2_4,
                colnames(result_MBCdeg2_5), result_MBCdeg2_5)
write.table(tmp, out3, sep="\t", append=F, quote=F, row.names=F, col.names=F)
