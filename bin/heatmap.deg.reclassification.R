#script by Raül García Veiga

#heatmap not classified
setwd("/no_backup/rg/ragarcia/repeated")

nclassHrt <- read.table(file='Hrt/expression/excluded/selected.genes.not.class.tsv',sep='\t', check.names=FALSE)
row.names(nclassHrt)=nclassHrt$V1
nclassHrt$V1 <- NULL
metaHrt <- read.table(file='Hrt/expression/excluded/metadata.tsv', sep='\t', check.names=FALSE)
colnames(metaHrt) <- c('names',	'class',	'V3')
metaHrt <- metaHrt[order(metaHrt$class), ]

znclassHrt <- t(apply(nclassHrt, 1, scale))
znclassHrt <- as.data.frame(znclassHrt)

colnames(znclassHrt) <- c('10.5',	'11.5',	'12.5',	'13.5',	'14.5',	'15.5',	'16.5',	'PN')
znclassHrt$extreme <- colnames(znclassHrt)[apply(znclassHrt,1,which.max)]

znclassHrt <- znclassHrt[complete.cases(znclassHrt), ]
znclassHrt <- znclassHrt[order(znclassHrt$extreme), ]

znclassHrt$names <- row.names(znclassHrt)
metaHrt <- merge(x = znclassHrt, y = metaHrt, by = "names")
znclassHrt$names <- NULL
znclassHrt$extreme <- NULL
row.names(metaHrt)=metaHrt$names
metaHrt <- select(metaHrt, 'class', 'extreme')
metaHrt <- metaHrt[order(metaHrt$extreme), ]
metaHrt$class <- NULL
znclassHrt.clust <- pheatmap(znclassHrt, main = 'Hrt not classified genes',cluster_cols = FALSE,  clustering_method = 'ward.D2', show_colnames = T, show_rownames = F, annotation_row = metaHrt)
znclassHrt.clust <- cbind(cluster = cutree(znclassHrt.clust$tree_row, k =  6))

clust.frame <- as.data.frame(znclassHrt.clust)
clust.frame$cluster <- as.factor(clust.frame$cluster)
pheatmap(znclassHrt, main = 'Hrt not classified genes', cluster_cols = FALSE, treeheight_row = 0,  cutree_rows = 6, clustering_method = 'ward.D2', show_colnames = T, show_rownames = F, annotation_row = clust.frame[,'cluster', drop = F])

matches <- c('bending','downregulation','downregulation','peaking', 'upregulation', 'peaking') #naming following cluster names for assigning to the output matrix
names(matches) <- c(1,2,3,4,5,6)

clust.frame$class <- matches[clust.frame$cluster]
write.table(clust.frame, "Hrt/expression/excluded/clusterinfo.tsv", sep="\t", quote = F)
