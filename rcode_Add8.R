setwd("‾/2020/研究/長部さん/")
rm(list = ls())
in_f <- "Additional8_pre.txt"

data <- read.table(in_f, header=TRUE, row.names=1, sep="¥t", quote="")

################
###  一致度：パターンは何でもよいので一致していればよい
###  table_002
################
hoge <- data
common <- matrix(10, nrow=ncol(hoge), ncol=ncol(hoge))
for(i in 1:ncol(hoge)){
  for(j in 1:ncol(hoge)){
    common[i,j] <- sum(hoge[,i] == hoge[,j])/nrow(hoge)
  }
}
namae <- colnames(data)
colnames(common) <- namae
rownames(common) <- namae
common
tmp <- cbind(namae, common)
write.table(tmp, "Additional8.txt", sep="¥t", append=F, quote=F, row.names=F)
