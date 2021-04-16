#boxplot violinplot dotplot for correlations 


rm(list = ls())

getwd()
setwd("/no_backup/rg/ragarcia/MBrain/H3K4me3")

library(ggplot2)


expression <- read.table('../expression/MBrain/expression.tsv', sep="\t" , header=T, row.names = NULL)
chrom <- read.table('H3K4me3.matrix.after.QN.merged.upreg.tsv', sep="\t")

expression <- subset(expression, row.names %in% chrom$V1)

row.names(chrom)=chrom$V1
chrom$V1 <- NULL

row.names(expression)=expression$row.names
expression$row.names <- NULL

#row.names(expression)=expression$gene_id
#expression$gene_id <- NULL

#row.names(expression)=expression$V1
#expression$V1 <- NULL


#ordering by rownames
chrom2 <- chrom[rownames(expression), ]

correl <- sapply(1:ncol(expression), function(i) cor(expression[,i], chrom2[,i], method = 'pearson'))
correl <- as.data.frame(correl)
time <- c('10','11','12','13','14','15','16','17')
correl <- cbind(correl, time)


texpression <- t(expression)
tchrom2 <- t(chrom2)
gencorrel <- sapply(1:ncol(texpression), function(i) cor(texpression[,i], tchrom2[,i], method = 'pearson'))

gencorrel <- as.data.frame(gencorrel)
medianWithoutNA<-function(x) {
  median(x[which(!is.na(x))])
}

correl$time <- as.numeric(correl$time)
class(correl$time)

ggplot() +
  geom_violin(data = gencorrel, 
              aes(x=4.5, y=gencorrel), 
              alpha=.3, colour="white", width = 8) +
  geom_boxplot(data = gencorrel, 
               aes(x=4.5, y=gencorrel),
               alpha = .4, outlier.shape = NA, width = 3) +
  geom_point(data = correl,
             aes(x=time, y=correl, size=time),
             shape = 21) +
  ylab(expression(paste("Correlation ", italic("(Pearson)")))) + xlab('H3K4me3') +
  scale_size_continuous(breaks = c(10, 11, 12, 13, 14, 15, 16, 17)) + labs(title = "Upreg > 25 TPM MBrain")

median(correl$correl) # median steady state
apply(gencorrel, 2, medianWithoutNA) #median time course












###########
###OLD CODE
###########
#chrom2
meta <- as.data.frame(rownames(chrom2))

rownames(meta) <- meta$`rownames(chrom2)`
meta$extreme <- colnames(chrom2)[apply(chrom2,1,which.min)]
meta$`rownames(chrom2)` <- NULL

zchrom2 <- t(apply(chrom2, 1, scale))
zchrom2 <- as.data.frame(zchrom2)

#colnames(zupMBrain) <- c('10.5',	'11.5',	'12.5',	'13.5',	'14.5',	'15.5',	'16.5',	'PN')
#zchrom2$extreme <- colnames(zchrom2)[apply(zchrom2,1,which.min)]

zchrom2 <- zchrom2[complete.cases(zchrom2), ]
zchrom2 <- zchrom2[order(zchrom2$extreme), ]

metaMBrain$names <- row.names(metaMBrain)
zupMBrain$names <- row.names(zupMBrain)
metaMBrain <- merge(x = zupMBrain, y = metaMBrain, by = "names")
zupMBrain$names <- NULL
zupMBrain$extreme <- NULL
row.names(metaMBrain)=metaMBrain$names
metaMBrain <- select(metaMBrain, 'class', 'extreme')
metaMBrain <- metaMBrain[order(metaMBrain$extreme), ]
metaMBrain$class <- NULL
pheatmap(zupMBrain, main = 'MBrain Upregulated genes',cluster_cols = FALSE, cluster_rows=F, clustering_method = 'ward.D2', show_colnames = T, show_rownames = F, annotation_row = metaMBrain)
