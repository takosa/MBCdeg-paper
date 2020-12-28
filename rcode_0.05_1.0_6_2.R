param_trial <- 50                      # number of trials
param_clust_num <- 2                   # number of clusters

param_G1 <- 6                          # number of samples in G1 group
param_G2 <- 6                          # number of samples in G2 group
param_Ngene <- 10000                   # number of genes
param_PDEG <- 0.05                     # fraction of DEGs
param_FC <- 4                          # fold change of DEGs
param_PG1 <- 1.0                       # fraction of PG1 in DEGs

library(TCC)
library(MBCluster.Seq)
library(ROC)

matome1 <- NULL
matome2 <- NULL
matome3 <- NULL
for(i in 1:param_trial){
print(i)
##########################################
###  Generation of simulation data
###  http://www.iu.a.u-tokyo.ac.jp/‾kadota/r_seq.html#count_simulation_RNAseq_biological_2_kiso_tcc
##########################################

set.seed(i)
tcc <- simulateReadCounts(Ngene=param_Ngene,
             PDEG=param_PDEG,
             DEG.assign=c(param_PG1, 1-param_PG1),
             DEG.foldchange=c(param_FC, param_FC),
             replicates=c(param_G1, param_G2))
data <- tcc$count
data.cl <- c(rep(1, param_G1), rep(2, param_G2))

tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=0.05)
result <- getResult(tcc, sort=FALSE)

ranking <- result$rank
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))

norm.factors <- tcc$norm.factors
ef.libsizes <- colSums(data)*norm.factors
size.factors <- ef.libsizes/mean(ef.libsizes)

auc_tcc <- auc

##########################################
###  MBCluster.Seq with default normalization
###  http://www.iu.a.u-tokyo.ac.jp/‾kadota/r_seq.html#analysis_clustering_RNAseq_genes_kiso_MBCluster.Seq
##########################################
hoge <- RNASeq.Data(data, Normalizer=NULL,
           Treatment=data.cl, GeneID=rownames(data))

c0 <- KmeansPlus.RNASeq(data=hoge, nK=param_clust_num,
           model="nbinom", print.steps=F)
#print(c("c0:", sprintf("%1.4f", tmp)), quote=F)

cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)

print(c("L2norm:", sprintf("%1.4f", L2norm)), quote=F)
print(res_table)

ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))


L2norm_def <- L2norm
auc_mbc_def <- auc
res_table_def <- res_table
k_def <- k

##########################################
###  MBCluster.Seq with TCC normalization
###  http://www.iu.a.u-tokyo.ac.jp/‾kadota/r_seq.html#analysis_clustering_RNAseq_genes_advanced_TCC_MBCluster.Seq
##########################################
hoge <- RNASeq.Data(data, Normalizer=log2(size.factors),
           Treatment=data.cl, GeneID=rownames(data))


c0 <- KmeansPlus.RNASeq(data=hoge, nK=param_clust_num,
           model="nbinom", print.steps=F)
#print(c("c0:", sprintf("%1.4f", tmp)), quote=F)


cls <- Cluster.RNASeq(data=hoge, model="nbinom",
           centers=c0$centers, method="EM")
res_table <- table(cls$cluster)
L2norm <- sqrt(rowSums(abs(cls$centers)^2))
k <- which.min(L2norm)

print(c("L2norm:", sprintf("%1.4f", L2norm)), quote=F)
print(res_table)


ranking <- rank(cls$probability[,k])
obj <- as.numeric(tcc$simulation$trueDEG != 0)
auc <- AUC(rocdemo.sca(truth=obj, data=-ranking))


L2norm_tcc <- L2norm
auc_mbc_tcc <- auc
res_table_tcc <- res_table
k_tcc <- k


matome1 <- rbind(matome1, c(auc_tcc, auc_mbc_def, auc_mbc_tcc))
matome2 <- rbind(matome2, c(L2norm_def, L2norm_def[k_def], L2norm_tcc, L2norm_tcc[k_tcc]))
matome3 <- rbind(matome3, c(res_table_def, res_table_tcc))
}


colnames(matome1) <- c("TCC", "MBC_def", "MBC_tcc")
matome1
out_f1 <- paste("res_",param_PDEG, "_", param_PG1, "_", param_G1, "_", param_clust_num, "_AUC.txt", sep="")
write.table(matome1, out_f1, sep="¥t", append=F, quote=F, row.names=F)
summary(matome1)


colnames(matome2) <- c("1_1", "1_2", "1_k", "2_1", "2_2", "2_k")
out_f2 <- paste("res_",param_PDEG, "_", param_PG1, "_", param_G1, "_", param_clust_num, "_L2norm.txt", sep="")
write.table(matome2, out_f2, sep="¥t", append=F, quote=F, row.names=F)
summary(matome2)


colnames(matome3) <- c("1_1", "1_2", "2_1", "2_2")
out_f3 <- paste("res_",param_PDEG, "_", param_PG1, "_", param_G1, "_", param_clust_num, "_cluster.txt", sep="")
write.table(matome3, out_f3, sep="¥t", append=F, quote=F, row.names=F)
summary(matome3)

