#########################################
### This file can be used to obtain the 
### raw AUC values showin in Figure 2 and 
### Additional file 2.
#########################################
###  Parameters  ###
N_trial <- 100             #number of trials
G <- 10000                 #number of genes
n1 <- 3                    #number of replicates for group 1
n2 <- 3                    #number of replicates for group 2
PDEG <- 0.05               #proportion of DEG
P1 <- 0.5                  #proportion of up-regulated DEGs in group 1
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
                 DEG.assign=c(P1, 1 - P1),
                 replicates=c(n1, n2))
tcc <- simulateReadCounts(Ngene=G,
             PDEG=PDEG,
             DEG.assign=c(P1, 1-P1),
             replicates=c(n1, n2),
             #fc.matrix=fc.matrix)
             DEG.foldchange=c(FC, FC))
data <- tcc$count
data.cl <- tcc$group$group
counts <- data
group <- data.cl

###  TCC  ###
#Obtaining size factors from TCC's normalization factors
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
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

#concatenate the results
matome <- rbind(matome, c(auc_MBCdeg1_2, auc_MBCdeg1_3, auc_MBCdeg1_4,
                          auc_MBCdeg2_2, auc_MBCdeg2_3, auc_MBCdeg2_4))
}

#Output(AUC)
colnames(matome) <- c("MBCdeg1(K=2)", "MBCdeg1(K=3)", "MBCdeg1(K=4)",
                      "MBCdeg2(K=2)", "MBCdeg2(K=3)", "MBCdeg2(K=4)")
out_f <- paste("Fig2_",PDEG, "_", P1, "_", n1, "_fixed.txt", sep="")
write.table(matome, out_f, sep="\t", append=F, quote=F, row.names=F)
#summary(matome)

