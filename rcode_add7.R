#########################################
### This file can be used to obtain the 
### raw AUC values showin in Additional file 7.
#########################################
###  Parameters  ###
N_trial <- 50              #number of trials
G <- 10000                 #number of genes
n1 <- 3                    #number of replicates for group 1
n2 <- 3                    #number of replicates for group 2
n3 <- 3                    #number of replicates for group 3
PDEG <- 0.25               #proportion of DEG
P1 <- 1/3                  #proportion of up-regulated DEGs in group 1
P2 <- 1/3                  #proportion of up-regulated DEGs in group 2
P3 <- 1/3                  #proportion of up-regulated DEGs in group 3
FC <- 4                    #degree of fold change (FC)

###  Load packages  ###
library(edgeR)
library(DESeq2)
library(TCC)
library(MBCluster.Seq)
library(ROC)

###################
###  Main loop  ###
###################
matome <- NULL
for(i in 1:N_trial){
print(i)

###  Generation of simulated data  ###
set.seed(i)
fc.matrix <- makeFCMatrix(Ngene=G,
                 PDEG=PDEG,
                 DEG.assign=c(P1, P2, P3),
                 replicates=c(n1, n2, n3))
tcc <- simulateReadCounts(Ngene=G,
             PDEG=PDEG,
             DEG.assign=c(P1, P2, P3),
             replicates=c(n1, n2, n3),
             fc.matrix=fc.matrix)          #Additional file 7
             #DEG.foldchange=c(FC, FC, FC))#Fig. 4 and Add. 6
data <- tcc$count
data.cl <- tcc$group$group
counts <- data
group <- data.cl

###  edgeR  ###
y <- edgeR::DGEList(counts = counts, group = group)
y <- edgeR::calcNormFactors(y)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
fit <- edgeR::glmQLFit(y, design)
coef <- 2:length(unique(data.cl))
qlf <- edgeR::glmQLFTest(fit, coef = coef)
res <- edgeR::topTags(qlf, n = nrow(counts), sort.by = "none")

#Calculation of AUC value
p.value <- res$table$PValue
ranking <- rank(p.value)
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_edger <- auc

###  DESeq2  ###
colData <- data.frame(condition=as.factor(data.cl))
d <- DESeqDataSetFromMatrix(countData=data,
             colData=colData, design=~condition)
d <- DESeq(d, test="LRT", full= ~condition, reduced= ~1)
tmp <- results(d)

#Calculation of AUC value
p.value <- tmp$pvalue
p.value[is.na(p.value)] <- 1
ranking <- rank(p.value)
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_deseq2 <- auc

###  TCC  ###
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=0.05)
result <- getResult(tcc, sort=FALSE)

#Calculation of AUC value
ranking <- result$rank
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_tcc <- auc

#Obtaining size factors from TCC's normalization factors
norm.factors <- tcc$norm.factors
ef.libsizes <- colSums(data)*norm.factors
size.factors <- ef.libsizes/mean(ef.libsizes)

###  MBCdeg1 (preparation)  ###
hoge <- RNASeq.Data(data, Normalizer=NULL,
           Treatment=data.cl, GeneID=rownames(data))

###  MBCdeg1 (K=2)  ###
K <- 2                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg1_2 <- auc

###  MBCdeg1 (K=3)  ###
K <- 3                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg1_3 <- auc

###  MBCdeg1 (K=4)  ###
K <- 4                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg1_4 <- auc

###  MBCdeg1 (K=5)  ###
K <- 5                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg1_5 <- auc

###  MBCdeg2 (preparation)  ###
hoge <- RNASeq.Data(data, Normalizer=log2(size.factors),
           Treatment=data.cl, GeneID=rownames(data))

###  MBCdeg2 (K=2)  ###
K <- 2                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg2_2 <- auc

###  MBCdeg2 (K=3)  ###
K <- 3                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg2_3 <- auc

###  MBCdeg2 (K=4)  ###
K <- 4                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg2_4 <- auc

###  MBCdeg2 (K=5)  ###
K <- 4                     #preselected number of clusters
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)            #number of genes assigned to each cluster
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster
#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc_MBCdeg2_5 <- auc

#concatenate the results
matome <- rbind(matome, c(auc_edger, auc_deseq2, auc_tcc,
                  auc_MBCdeg1_2, auc_MBCdeg1_3, auc_MBCdeg1_4, auc_MBCdeg1_5,
                  auc_MBCdeg2_2, auc_MBCdeg2_3, auc_MBCdeg2_4, auc_MBCdeg2_5))
}

#Output(AUC)
colnames(matome) <- c("edgeR", "DESeq2", "TCC",
                    "MBCdeg1(K=2)", "MBCdeg1(K=3)", "MBCdeg1(K=4)", "MBCdeg1(K=5)",
                    "MBCdeg2(K=2)", "MBCdeg2(K=3)", "MBCdeg2(K=4)", "MBCdeg2(K=5)")
matome
out_f <- paste("Add7_",PDEG, "_", sprintf("%3.2f", P1), "_", n1, "_gamma.txt", sep="")
write.table(matome, out_f, sep="\t", append=F, quote=F, row.names=F)
summary(matome)



