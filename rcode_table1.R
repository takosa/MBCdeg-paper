#########################################
### This file can be used to obtain the 
### contents showin in Table 1.
### If we use "set.seed(14)" instead of "set.seed(75)", 
### we will get similar results to Table 1b 
#########################################
###  Parameters  ###
G <- 10000                 #number of genes
n1 <- 3                    #number of replicates for group 1
n2 <- 3                    #number of replicates for group 2
PDEG <- 0.25               #proportion of DEG
P1 <- 0.5                  #proportion of up-regulated DEGs in group 1
FC <- 4                    #degree of fold change (FC)
K <- 3                     #preselected number of clusters

###  Load packages  ###
library(TCC)
library(MBCluster.Seq)
library(ROC)

##################
###  Table 1a  ###
##################
###  Generation of simulated data  ###
set.seed(1)
tcc <- simulateReadCounts(Ngene=G,
             PDEG=PDEG,
             DEG.assign=c(P1, 1-P1),
             replicates=c(n1, n2),
             DEG.foldchange=c(FC, FC))
data <- tcc$count
data.cl <- c(rep(1, n1), rep(2, n2))

###  TCC  ###
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
#Obtaining size factors from TCC's normalization factors
norm.factors <- tcc$norm.factors
ef.libsizes <- colSums(data)*norm.factors
size.factors <- ef.libsizes/mean(ef.libsizes)

###  MBCdeg2  ###
hoge <- RNASeq.Data(data, Normalizer=log2(size.factors),
           Treatment=data.cl, GeneID=rownames(data))
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster

#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))

#concatenate the results
out <- matrix(0, nrow=3, ncol=5)
hoge <- cls$cluster[1:(G*PDEG*P1)]         # DEG1
for(i in 1:K){out[i,1] <- sum(hoge == i)}  # DEG1
if(P1 == 1.0){hoge <- NULL                # DEG2
}else{hoge <- cls$cluster[(G*PDEG*P1+1):(G*PDEG)]}# DEG2
for(i in 1:K){out[i,2] <- sum(hoge == i)}  # DEG2
hoge <- cls$cluster[(G*PDEG+1):G]          # non-DEG
for(i in 1:K){out[i,3] <- sum(hoge == i)}  # non-DEG
hoge <- cls$cluster                        # Total(all genes)
for(i in 1:K){out[i,4] <- sum(hoge == i)}  # Total(all genes)
out[,5] <- L2norm
colnames(out) <- c("DEG1", "DEG2", "nonDEG", "Total", "L2Norm")
Table1a <- out
AUC_Table1a <- auc

##################
###  Table 1b  ###
##################
###  Generation of simulated data  ###
set.seed(75)
tcc <- simulateReadCounts(Ngene=G,
             PDEG=PDEG,
             DEG.assign=c(P1, 1-P1),
             replicates=c(n1, n2),
             DEG.foldchange=c(FC, FC))
data <- tcc$count
data.cl <- c(rep(1, n1), rep(2, n2))

###  TCC  ###
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
#Obtaining size factors from TCC's normalization factors
norm.factors <- tcc$norm.factors
ef.libsizes <- colSums(data)*norm.factors
size.factors <- ef.libsizes/mean(ef.libsizes)

###  MBCdeg2  ###
hoge <- RNASeq.Data(data, Normalizer=log2(size.factors),
           Treatment=data.cl, GeneID=rownames(data))
c0 <- KmeansPlus.RNASeq(data=hoge, nK=K,
           model="nbinom", print.steps=F)
cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
L2norm <- sqrt(rowSums(abs(cls$centers)^2))#calculation of L2 Norm
k <- which.min(L2norm)                     #index for the non-DEG cluster

#Calculation of AUC value
ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))
auc

#concatenate the results
out <- matrix(0, nrow=3, ncol=5)
hoge <- cls$cluster[1:(G*PDEG*P1)]         # DEG1
for(i in 1:K){out[i,1] <- sum(hoge == i)}  # DEG1
if(P1 == 1.0){hoge <- NULL                # DEG2
}else{hoge <- cls$cluster[(G*PDEG*P1+1):(G*PDEG)]}# DEG2
for(i in 1:K){out[i,2] <- sum(hoge == i)}  # DEG2
hoge <- cls$cluster[(G*PDEG+1):G]          # non-DEG
for(i in 1:K){out[i,3] <- sum(hoge == i)}  # non-DEG
hoge <- cls$cluster                        # Total(all genes)
for(i in 1:K){out[i,4] <- sum(hoge == i)}  # Total(all genes)
out[,5] <- L2norm
colnames(out) <- c("DEG1", "DEG2", "nonDEG", "Total", "L2Norm")
Table1b <- out
AUC_Table1b <- auc

##################
###  Output    ###
##################
AUC_Table1a
Table1a
AUC_Table1b
Table1b
