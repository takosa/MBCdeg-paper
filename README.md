# MBCdeg-paper

This repository contains some code used in the paper "Differential expression analysis using a model-based gene clustering algorithm for RNA-seq data(under review)"

This repository contains a total of 16 files: 12 R-code files (*.R), one sample data file (sample.txt), the R-code files for executing the sample file (sample_MBCdeg1.R and sample_MBCdeg2.R), and me. Followings are the details for individual files. You won't get exactly the same results, but you will get similar results. 

![](method_description.svg)

###  MBCdeg.R  ###
This file contains MBCdeg function which execute MBCdeg procedure. This function dependent `MBCluster.Seq` packages.

Usage:

```r
source("MBCdeg.R")
count_matrix <- read.table("sample.txt")
treatment <- data.cl <- c(rep(1, 5), rep(2, 6))
no_clusters <- 3
result <- MBCdeg(counts = count_matrix, treatment = treatment, K = no_clusters)
print(result)
```

###  rcode_fig1.R  ###
By executing this file with the default parameter settings (i.e., 100 trials, G = 10,000, n1 = n2 = 3, PDEG = 0.05, P1 = 0.5, FC = 4, and K = 3), one can obtain a tab-delimited file (named "Fig1_0.05_0.5_3_fixed.txt") that contains raw AUC values of five methods for individual trials under PDEG = 0.05 and P1 = 0.5 in Figure 1. One can obtain all the raw AUC values in this figure by changing the two parameters (i.e., PDEG and P1). One can also obtain the raw AUC values in Additional file 1 by changing the four parameters (n1, n2, PDEG, and P1).

###  rcode_fig2.R  ###
By executing this file with the default parameter settings (i.e., 100 trials, G = 10,000, n1 = n2 = 3, PDEG = 0.05, P1 = 0.5, and FC = 4), one can obtain a tab-delimited file (named "Fig2_0.05_0.5_3_fixed.txt") that contains raw AUC values of MBCdeg with K=2-4 for individual trials under PDEG = 0.05 and P1 = 0.5 in Figure 2. One can obtain all the raw AUC values in this figure by changing the two parameters (i.e., PDEG and P1). One can also obtain the raw AUC values in Additional file 2 by changing the four parameters (n1, n2, PDEG, and P1).

###  rcode_table1.R  ###
By executing this file, one can obtain the contents shown in Table 1.

###  rcode_add3.R  ###
By executing this file with the default parameter settings (i.e., 100 trials, G = 10,000, n1 = n2 = 3, PDEG = 0.05, P1 = 0.5, FC = 4, and K = 3), one can obtain a tab-delimited file (named "Add3_0.05_0.5_3_gamma.txt") that contains raw AUC values of five methods for individual trials under n1 = n2 = 3, PDEG = 0.05, and P1 = 0.5 in Additional file 3. One can obtain all the raw AUC values in Additional file 3 by changing the four parameters (n1, n2, PDEG, and P1).

###  rcode_add4.R  ###
By executing this file with the default parameter settings (i.e., 100 trials, G = 10,000, n1 = n2 = 3, PDEG = 0.05, P1 = 0.5, and FC = 4), one can obtain a tab-delimited file (named "Add4_0.05_0.5_3_gamma.txt") that contains raw AUC values of MBCdeg with K=2-4 for individual trials under n1 = n2 = 3, PDEG = 0.05, and P1 = 0.5 in Additional file 4. One can obtain all the raw AUC values in Additional file 4 by changing the four parameters (n1, n2, PDEG, and P1).

###  rcode_fig3.R  ###
By executing this file with the default parameter settings (i.e., 50 trials, G = 10,000, n1 = n2 = 3, PDEG = 0.45, P1 = 0.5, FC = 4, and K = 3), one can obtain a tab-delimited file (named "Fig3_0.45_0.5_3_fixed.txt") that contains raw AUC values of five methods for individual trials under PDEG = 0.45 and P1 = 0.5 in Figure 3. One can obtain all the raw AUC values in this figure by changing the two parameters (i.e., PDEG and P1). One can also obtain the raw AUC values in Additional file 5 by changing the four parameters (n1, n2, PDEG, and P1).

###  rcode_table2.R  ###
By executing this file, one can obtain the contents shown in Table 2.

###  rcode_table3.R  ###
By executing this file, one can obtain the contents shown in Table 3.

###  rcode_fig4.R  ###
By executing this file with the default parameter settings (i.e., 50 trials, G = 10,000, n1 = n2 = n3 = 3, PDEG = 0.25, P1 = P2 = P3 = 1/3, and FC = 4), one can obtain a tab-delimited file (named "Fig4_0.25_0.33_3_fixed.txt") that contains raw AUC values of five methods for individual trials under (1/3, 1/3, 1/3) in Figure 4. One can obtain all the raw AUC values in this figure by changing the three parameters (i.e., P1, P2, and P3). One can also obtain the raw AUC values in Additional file 6 by changing a total of six parameters (n1, n2, n3, P1, P2, and P3).

###  rcode_add7.R  ###
By executing this file with the default parameter settings (i.e., 50 trials, G = 10,000, n1 = n2 = n3 = 3, PDEG = 0.25, P1 = P2 = P3 = 1/3, and FC = 4), one can obtain a tab-delimited file (named "Add7_0.25_0.33_3_gamma.txt") that contains raw AUC values of five methods for individual trials under n1 = n2 = n3 = 3 and (1/3, 1/3, 1/3) in Additional file 7. One can obtain all the raw AUC values in this figure by changing a total of six parameters (n1, n2, n3, P1, P2, and P3).

###  rcode_table4.R  ###
By executing this file, one can obtain the contents shown in Table 4 and Additional file 8.

###  sample.txt  ###
This file is a template for the input file that the user gives when actually using MBCdeg. This file was generated by the simulateReadCounts function in TCC with following parameters: G = 2000, n1 = 5, n2 = 6, PDEG = 0.2, P1 = 0.9 (i.e., P2 = 0.1), FC1 = 4 for up-regulated in group 1, and FC2 = 9 for up-regulated in group 2. Accordingly, 2000*0.2 = 400 genes were designed as DEGs. Of these, the first 90% (i.e., 400*0.9 = 360 DEG1 genes) were up-regulated in group 1 and the remaining 10% (400*0.1 = 40 DEG2 genes) were up-regulated in group 2. Although the latter gene number is small, the degree of FC (FC2 = 9) is much larger than that of the former (i.e., FC1 = 4); We hypothesized that the 40 DEG2 genes could also form distinct clusters. Therefore the truth of this dataset is as follows: the first 360 genes (gene_1, gene_2, ..., and gene_360) have the DEG1 pattern, the next 40 genes (gene_361, gene_362, ..., and gene_400) have the DEG2 pattern, and the remaining 1,600 genes (gene_401, gene_402, ..., gene_2000) have the non-DEG pattern.

###  sample_MBCdeg1.R  ###
This is an exammple R code to analyze the two-group sample data ("sample.txt") by using MBCdeg1 with K = 3. The input data compares group A (5 replicates) vs. group B (6 replicates). By executing this file with the default parameter settings, one can obtain a tab-delimited text file (named "sample_MBCdeg1.txt") as the output. The output file has seven columns of information added to the right side of the input file. The first three columns (named "1", "2", and "3") correspond to the posterior probabilities (PPs) assigned to individual clusters. The next column (named "PP_nonDEG") corresponds to the PPs determined as non-DEG cluster. In this case, we see the third cluster (named "3") was identified as the non-DEG cluster. The fifth column (named "ranking") corresponds to overall gene ranking that is obtained based on the non-DEG PPs. The sixth column (named "pat") corresponds to the pattern name assigned to each gene (either "DEG1", "DEG2" or "non-DEG"). The information of the last column (named "clust_num") corresponds to the cluster name (either "1", "2" or "3")is essentially the same as that of the sixth column. In this case, the first, second, and third clusters correspond to DEG2, DEG1, and non-DEG patterns, respectively. The information displayed on the "R Console" screen may also be useful.

###  sample_MBCdeg2.R  ###
This is an exammple R code to analyze the two-group sample data ("sample.txt") by using MBCdeg2 with K = 3. The others are essentially the same as those described above.



