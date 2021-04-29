setwd("/no_backup/rg/ragarcia/MBrain/H3K9ac")

library(ggplot2)

expression <- read.table('../expression/MBrain/expression.r10.tsv', sep="\t" , header=T, row.names = NULL)
chrom <- read.table('H3K9ac.matrix.after.QN.merged.tsv', sep="\t" )
#expression <- subset(expression, row.names %in% chrom$V1)

#row.names(chrom)=chrom$V1
#chrom$V1 <- NULL

row.names(expression)=expression$row.names
expression$row.names <- NULL

#row.names(expression)=expression$gene_id
#expression$gene_id <- NULL

#ordering by rownames
chrom2 <- chrom[rownames(expression), ]

correl <- sapply(1:ncol(expression), function(i) cor(expression[,i], chrom2[,i], method = 'pearson'))
correl <- as.data.frame(correl)
time <- c('10','11','12','13','14','15','16')
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
  ylab(expression(paste("Correlation ", italic("(Pearson)")))) + xlab('H3K9ac') +
  scale_size_continuous(breaks = c(10, 11, 12, 13, 14, 15, 16)) + labs(title = 'All genes MBrain')

median(correl$correl) # median time course
apply(gencorrel, 2, medianWithoutNA) #median steady state
