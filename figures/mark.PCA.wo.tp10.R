setwd("/no_backup/rg/ragarcia/MBrain/H3K9ac")
library(ggplot2)

chrom <- read.table('H3K9ac.matrix.after.QN.merged.tsv', sep="\t")
row.names(chrom)=chrom$V1
chrom$V1 <- NULL
chrom.subset <- read.table('H3K9ac.matrix.after.QN.merged.tsv', sep="\t")

colnames(chrom) <- c('X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN')
pca_res <- prcomp(t(chrom), scale = F)
summary(pca_res)
df <- as.data.frame(pca_res$x)
df$tp <- rownames(df)

ggplot(df, aes(x=PC1, y=PC2, size=as.factor(tp))) +
  geom_point() + ggtitle("H3K9ac MBrain") + xlab(summary(pca_res)$importance[2,1]*100) +
  ylab(summary(pca_res)$importance[2,2]*100) +
  geom_label(
    label=colnames(chrom), 
    nudge_x = 20) + 
  theme(legend.position = "none")
