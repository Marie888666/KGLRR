# KGLRR_
# KGLRR: A low-rank representation K-means with graph regularization constraint method for Single-cell type identification given 


KGLRR  Overview:



KGLRR performs simultaneously dimensionality reduction and clustering using low-rank representation and K-means clustering. 

KGLRR is written in the MATLAB programming language. To use, please download the KGLRR folder and follow the instructions provided in the “README.doc”.

Files:
KGLRR.m---A script with a real scRNA-seq data to shows how to run the code.

Data
data$Treutlin- Treutlin dataset original expression data. The original data can be downloaded from GEO website. Login number is GSE52583

Preprocessing approach
FilterGenesZero.m---The gene filter step, we remove the genes whose expressions (the expression of the gene is nonzero) are <5% of all cells.

normalize.m---L2-norm is applied on gene expression of each cell to eliminate the scale differences among samples.

KGLRR approach
lrr_relaxed.m---The core code of KGLRR. KGLRR performs simultaneously dimensionality reduction and clustering using low-rank representation and K-means clustering. 




Please send any questions or found bugs to wanglinpingqfu@126.com.
