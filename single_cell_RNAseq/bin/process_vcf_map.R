##
## process_vcf_map.R
## Adapted from SICCA pipeline code from Ravi Patel and Kim Taylor
##
## determine high-confidence matches and plot results
##
## input: matrix from vcf-match-sample-ids
## output: return 0 if all good, 1 if issues; leave plots in results directory

# below is visible in singularity env
library(pheatmap)


args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Please provide pool")
}
pool=args[1]

matrix <- as.matrix(read.table(mapFile, header=TRUE, row.names=1, sep ="\t"))

# z standardization by rows
matrix.zR <- t(scale(t(matrix),center=TRUE,scale=TRUE))
# 1 if z > 1.96 (p<0.05)
matrix.1R <- (matrix.zR > 1.96)*1

# z standardization by cols
matrix.zC <- scale(matrix,center=TRUE,scale=TRUE)
matrix.1C <- (matrix.zC > 1.96)*1

# if single outlier by column and sample agree, succeed
# and the #1s = #samples (row dimension)
result = identical(matrix.1C,matrix.1R) & (sum(matrix.1C) == dim(matrix)[1])

## leave plots regardless
pdf(paste(mapDir,"VCFmap_heatmap.pdf",sep=""))
title=sprintf("Pool %s pairwise likelihoods.",i)
pheatmap(matrix, clustering_method = "ward.D2", cluster_cols = T, main = title)
title=sprintf("Pool %s column-wise z-scores.",i)
pheatmap(matrix.zC, clustering_method = "ward.D2", cluster_cols = T, main = title)
title=sprintf("Pool %s row-wise z-scores.",i)
pheatmap(matrix.zR, clustering_method = "ward.D2", cluster_cols = T, main = title)
dev.off()

if (result == FALSE) {
	print("process_vcf_map.R: mapping not confidently determined.")
	print("Z scores by column:")
	print(matrix.zC)
	print("Z scores by row:")
	print(matrix.zR)
} 

## output the map regardless
# I'm sure there's a more elegant way to do this...
map<-as.data.frame(matrix(nrow=dim(matrix)[1],ncol=2))
colnames(map)<-c("SampleID","FMcluster")
map$SampleID = rownames(matrix.1C)
for (i in 1:dim(matrix)[1]) {
	for (j in 1:dim(matrix)[2]) {
		if (matrix.1C[i,j] == 1 & matrix.1R[i,j] == 1) {
			map[i,"FMcluster"] = colnames(matrix.1C)[j]
		}
	}
}
write.csv(map,file=paste(mapDir,"VCFmap.csv",sep=""),row.names=FALSE,quote=FALSE)


quit(status=1-result)




