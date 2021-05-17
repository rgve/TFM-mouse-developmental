expr <- read.table(file='/no_backup/rg/ragarcia/repeated/MBrain/expression/excluded/selected.genes.expr.0.upreg.tsv',sep='\t', check.names=FALSE)
row.names(expr)=expr$V1
expr$V1 <- NULL
metaHrt <- read.table(file='/no_backup/rg/ragarcia/repeated/MBrain/expression/excluded/metadata.tsv', sep='\t', check.names=FALSE)
row.names(metaHrt)=metaHrt$V1
metaHrt$V1 <- NULL
colnames(metaHrt) <- c('class',	'time_point')
metaHrt <- metaHrt[order(metaHrt$class), ]

expr <- t(apply(expr, 1, scale))
colnames(expr) <- c('10.5',	'11.5',	'12.5',	'13.5',	'14.5',	'15.5',	'16.5',	'PN')
expr <- as.data.frame(expr)

expr$extreme <- colnames(expr)[apply(expr,1,which.max)]

expr <- expr[complete.cases(expr), ]
expr <- expr[order(expr$extreme), ]

metaHrt$names <- row.names(metaHrt)
expr$names <- row.names(expr)
metaHrt <- merge(x = expr, y = metaHrt, by = "names")
expr$names <- NULL
expr$extreme <- NULL
row.names(metaHrt)=metaHrt$names
metaHrt <- select(metaHrt, 'class', 'extreme')
metaHrt <- metaHrt[order(metaHrt$extreme), ]
metaHrt$class <- NULL
pheatmap(expr, cluster_cols = FALSE, cluster_rows=F, main = 'Hrt Upregulated genes',clustering_method = 'ward.D2', show_colnames = T, show_rownames = F, annotation_row = metaHrt)
