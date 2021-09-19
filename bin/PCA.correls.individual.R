#script by Raül García Veiga

library(ggplot2)
#for -10 tp just delete tp 10 tags


#merged mark matrix
chrom <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/results/marks/H3K27me3/H3K27me3.matrix.after.QN.merged.tsv', sep="\t")
colnames(chrom) <- c('10.5','11.5', '12.5', '13.5', '14.5', '15.5', '16.5', 'PN')

#R2 mark matrix
chrom2 <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/results/marks/H3K27me3/H3K27me3.R2.matrix.after.QN.merged.tsv', sep="\t")
colnames(chrom2) <- c('10.2','11.2', '12.2', '13.2', '14.2', '15.2', '16.2', 'PN.2')

#R1 mark matrix
chrom3 <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/results/marks/H3K27me3/H3K27me3.R1.matrix.after.QN.merged.tsv', sep="\t")
colnames(chrom3) <- c('10.1','11.1', '12.1', '13.1', '14.1', '15.1', '16.1', 'PN.1')

#diagonal correl for barplot
co <- diag(cor(chrom3,chrom2, method = 'spearman'))

#matrix correl for heatmap
corep <- cor(chrom3,chrom2, method = 'spearman')

co <- as.data.frame(co)
co <- t(co)
colnames(co) <- c('10', '11', '12', '13', '14', '15', '16', 'PN')

barplot((co), ylim = c(0,1) ,main = 'H3K27me3 rep. correlation', col='coral',xlab = 'spearman')

pheatmap(corep, main = 'H3K27me3 rep. r spearman', show_rownames = T, cluster_cols = FALSE, cluster_rows =F)

#PCA merged matrix
pca_res <- prcomp(t(chrom), scale = F)
summary(pca_res)
df <- as.data.frame(pca_res$x)
df$tp <- rownames(df)

ggplot(df, aes(x=PC1, y=PC2, size=as.factor(tp))) +
  geom_point() + ggtitle("H3K27me3 MBrain") + xlab(summary(pca_res)$importance[2,1]*100) +
  ylab(summary(pca_res)$importance[2,2]*100) +
  geom_label(
    label=colnames(chrom), 
    nudge_x = 20) + 
  theme(legend.position = "none")

#PCA r1 
pca_res1 <- prcomp(t(chrom3), scale = F)
summary(pca_res1)
df <- as.data.frame(pca_res1$x)
df$tp <- rownames(df)

ggplot(df, aes(x=PC1, y=PC2, size=as.factor(tp))) +
  geom_point() + ggtitle("H3K27me3 MBrain R1") + xlab(summary(pca_res1)$importance[2,1]*100) +
  ylab(summary(pca_res1)$importance[2,2]*100) +
  geom_label(
    label=colnames(chrom), 
    nudge_x = 20) + 
  theme(legend.position = "none")

#PCA r2
pca_res2 <- prcomp(t(chrom2), scale = F)
summary(pca_res2)
df <- as.data.frame(pca_res2$x)
df$tp <- rownames(df)

ggplot(df, aes(x=PC1, y=PC2, size=as.factor(tp))) +
  geom_point() + ggtitle("H3K27me3 MBrain R2") + xlab(summary(pca_res2)$importance[2,1]*100) +
  ylab(summary(pca_res2)$importance[2,2]*100) +
  geom_label(
    label=colnames(chrom), 
    nudge_x = 20) + 
  theme(legend.position = "none")
