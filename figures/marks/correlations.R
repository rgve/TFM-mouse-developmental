#boxplot violinplot dotplot and correlations
#for marks without first timepoint delete (2) time tags #17 and remove first expression matrix column
library(ggplot2)

expression <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/expression/excluded/selected.genes.expr.25.upreg.tsv', sep="\t" , header=F)
rownames(expression) <- expression$V1
expression$V1 <-  NULL

#delete after '.' in rownames for matching mark matrix
rownames(expression) <- sub("[.].*", "", rownames(expression))

chrom <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/results/marks/H3K36me3/H3K36me3.matrix.after.QN.merged.tsv', sep="\t")

#delete first expression timepoint for marks without tp 10.5
  #expression$X10.5.1 <- NULL
  #expression$V2 <- NULL


#get only chrom rows present in expression subset
colnames(expression) <- colnames(chrom)

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
  ylab(expression(paste("Correlation ", italic("(Pearson)")))) + xlab('H3K36me3') +
  scale_size_continuous(breaks = c(10, 11, 12, 13, 14, 15, 16, 17)) + labs(title = 'Upreg >25 M-Brain')
                    
# get time course median
median(correl$correl) 
                    
# get steady state median
apply(gencorrel, 2, medianWithoutNA)
