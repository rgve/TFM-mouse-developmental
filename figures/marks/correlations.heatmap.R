setwd("/no_backup/rg/ragarcia/MBrain/H3K36me3")
library(ggplot2)

expression <- read.table('../expression/MBrain/expression.tsv', sep="\t" , header=T, row.names = NULL)
chrom <- read.table('H3K36me3.matrix.after.QN.merged.upreg25.tsv', sep="\t")
expression <- subset(expression, row.names %in% chrom$V1)

row.names(chrom)=chrom$V1
chrom$V1 <- NULL

row.names(expression)=expression$row.names
expression$row.names <- NULL

#row.names(expression)=expression$gene_id
#expression$gene_id <- NULL

#row.names(expression)=expression$V1
#expression$V1 <- NULL

#matching by rownames
chrom2 <- chrom[rownames(expression), ]
colnames(chrom2) <- c('10.5',	'11.5',	'12.5',	'13.5',	'14.5',	'15.5',	'16.5',	'PN')

#Manual ordering
chrom2$extreme <- colnames(chrom2)[apply(chrom2,1,which.max)]
chrom2 <- chrom2[complete.cases(chrom2), ]
chrom2 <- chrom2[order(chrom2$extreme), ]
chrom2$extreme <- NULL

colnames(expression) <- c('10.5',	'11.5',	'12.5',	'13.5',	'14.5',	'15.5',	'16.5',	'PN')
metaex <- as.data.frame(rownames(expression))
rownames(metaex) <- metaex$`rownames(expression)`
metaex$max.expr <- colnames(expression)[apply(expression,1,which.max)]
metaex$`rownames(expression)` <- NULL

expression2 <- t(apply(expression, 1, scale))
expression2 <- as.data.frame(expression2)
colnames(expression2) <- c('10.5',	'11.5',	'12.5',	'13.5',	'14.5',	'15.5',	'16.5',	'PN')
expression2$extreme <- colnames(expression2)[apply(expression2,1,which.max)]
expression2 <- expression2[complete.cases(expression2), ]
expression2 <- expression2[order(expression2$extreme), ]

zchrom <- t(apply(chrom2, 1, scale))
zchrom <- as.data.frame(zchrom)
colnames(zchrom) <- c('10.5',	'11.5',	'12.5',	'13.5',	'14.5',	'15.5',	'16.5',	'PN')
zchrom <- zchrom[rownames(expression2), ]

pheatmap(zchrom, cluster_cols = FALSE, cluster_rows=F, main = 'H3K36me3 >25 TPM',clustering_method = 'ward.D2', show_colnames = T, show_rownames = F, annotation_row = metaex )
