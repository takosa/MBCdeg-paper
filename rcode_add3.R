#########################################
### This file can be used to obtain the 
### raw AUC values showin in Additional file 3.
#########################################
####  Parameters  ###
N_trial <- 100             #number of trials
G <- 10000                 #number of genes
n1 <- 3                    #number of replicates for group 1
n2 <- 3                    #number of replicates for group 2
PDEG <- 0.05               #proportion of DEG
P1 <- 0.5                  #proportion of up-regulated DEGs in group 1
FC <- 4                    #degree of fold change (FC)
K <- 3                     #preselected number of clusters

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
             fc.matrix=fc.matrix)
             #DEG.foldchange=c(FC, FC))
data <- tcc$count
data.cl <- c(rep(1, n1), rep(2, n2))
counts <- data
group <- data.cl

###  edgeR  ###
y <- edgeR::DGEList(counts = counts, group = group)
y <- edgeR::calcNormFactors(y)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
fit <- edgeR::glmQLFit(y, design)
qlf <- edgeR::glmQLFTest(fit, coef = 2)
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
d <- DESeq(d)
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

###  MBCdeg1  ###
hoge <- RNASeq.Data(data, Normalizer=NULL,
           Treatment=data.cl, GeneID=rownames(data))
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

#rename
L2norm_MBCdeg1 <- L2norm
auc_MBCdeg1 <- auc
res_table_MBCdeg1 <- res_table
k_MBCdeg1 <- k

###  MBCdeg2  ###
hoge <- RNASeq.Data(data, Normalizer=log2(size.factors),
           Treatment=data.cl, GeneID=rownames(data))
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

#rename
L2norm_MBCdeg2 <- L2norm
auc_MBCdeg2 <- auc
res_table_MBCdeg2 <- res_table
k_MBCdeg2 <- k

#concatenate the results
matome <- rbind(matome, c(auc_edger, auc_deseq2, auc_tcc, auc_MBCdeg1, auc_MBCdeg2))
}

#Output(AUC)
colnames(matome) <- c("edgeR", "DESeq2", "TCC", "MBCdeg1", "MBCdeg2")
matome
out_f <- paste("Add3_",PDEG, "_", P1, "_", n1, "_gamma.txt", sep="")
write.table(matome, out_f, sep="\t", append=F, quote=F, row.names=F)
summary(matome)



