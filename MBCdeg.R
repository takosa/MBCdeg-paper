################################################################################
## This file contains functions to run MBCdeg procedure.
## Usage example:
##   source("MBCdeg.R")
##   count_matrix <- read.table("sample.txt")
##   treatment <- data.cl <- c(rep(1, 5), rep(2, 6))
##   no_clusters <- 3
##   result <- MBCdeg(counts = count_matrix, treatment = treatment, K = no_clusters)
##   print(result)
################################################################################

library(MBCluster.Seq)

#' Add together two numbers
#'
#' @param counts A count matrix
#' @param treatment A vector indicating information about experimental design
#' @param K A number of clusters
#' @param normalizer A vector which is passed to MBCluster.Seq::RNASeq.Data function
#' @param geneid A vector which is passed to MBCluster.Seq::RNASeq.Data function
#' @return A mbc object
MBCdeg <- function(counts, treatment, K, normalizer=NULL, geneid=NULL) {
  
  message("1. Preprocessing data...")
  mbc <- RNASeq.Data(counts, Normalizer=normalizer,
                     Treatment=treatment, GeneID=geneid)
  message("2. Initializing centers...")
  c0 <- KmeansPlus.RNASeq(data=mbc, nK=K, model="nbinom", print.steps=F)
  message("3. Clustering data...")
  capture.output({
    cls <- Cluster.RNASeq(data=mbc, model="nbinom", centers=c0$centers, method="EM")
  })
  
  ###  Identify the k value that corresponds to the non-DEG cluster  ###
  L2norm <- sqrt(rowSums(abs(cls$centers)^2))
  k <- which.min(L2norm)
  
  
  rtn <- list()
  rtn$centers <- cls$centers
  rtn$PP <- cls$probability
  rtn$PP_nonDEG <- cls$probability[,k]
  rtn$ranking <- rank(cls$probability[,k])
  rtn$cluster <- cls$cluster
  
  class(rtn) <- "mbc"
  rtn
}


##### Generic functions #####

#' print mbc object
#'
#' @param obj A mbc object
print.mbc <- function(obj) {
  cat("This is a mbc class object.\n")
  cat("It contains some values below.\n")
  cat("  centers: centers of each cluster.\n")
  cat("  PP: posterior probability.\n")
  cat("  PP_nonDEG: posterior probability of nonDEG cluster.\n")
  cat("  ranking: ranking of genes by PP_nonDEG.\n")
  cat("  cluster: cluster to which each gene is assigned.\n")
}
