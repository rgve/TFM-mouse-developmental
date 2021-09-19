#script by Raül García Veiga

#boxplot violinplot dotplot and correlations 

#setwd("/no_backup/rg/ragarcia/MBrain/H3K27ac")

library(ggplot2)


expression <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/expression/selected.genes.rep.1.2.log.after.QN.merged2.tsv', sep="\t" , header=T)#, row.names = NULL)

#expression <- read.table('/no_backup/rg/ragarcia/MBrain/expression/MBrain/selected.genes.rep.1.2.log.after.QN.merged2.tsv', sep="\t" , header=T)#, row.names = NULL)



chrom <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/results/marks/H3K27ac/H3K27ac.matrix.after.QN.merged.cut.tsv', sep="\t")


#chrom$V1 <- row.names(chrom)

#row.names(expression)=expression$row.n
#expression$V1 <- NULL

#row.names(expression)=expression$gene_id
#expression$gene_id <- NULL

#ordering by rownames
chrom2 <- chrom[rownames(expression), ]

colnames(expression) <- colnames(chrom2)


chrom2 <- as.matrix(chrom2)
expression <- as.matrix(expression)

#row.names(chrom)=chrom$V1
#chrom$V1 <- NULL

#row.names(expression)=expression$row.names
#expression$row.names <- NULL

#row.names(expression)=expression$gene_id
#expression$gene_id <- NULL

#row.names(expression)=expression$V1
#expression$V1 <- NULL
#row.names(expression)=expression$row.names
#expression$row.names <- NULL


#ordering by rownames
#chrom2 <- chrom[rownames(expression), ]

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
  ylab(expression(paste("Correlation ", italic("(Pearson)")))) + xlab('H3K27ac') +
  scale_size_continuous(breaks = c(10, 11, 12, 13, 14, 15, 16, 17)) + labs(title = "All genes MBrain")

median(correl$correl) # median steady state
apply(gencorrel, 2, medianWithoutNA) #median time course
