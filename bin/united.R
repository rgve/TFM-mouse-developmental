#!/usr/bin/Rscript

#script by Raül García Veiga

file_in <- commandArgs(trailingOnly = TRUE)

data_matrix <- read.table(file=file_in,sep='\t', check.names=FALSE, header=F)

data_matrix[, 2:17] <- log2(data_matrix[2:17]+1)

row.names(data_matrix)=data_matrix$V1
data_matrix$V1=NULL

colnames(data_matrix)=c('10.5.1',	'11.5.1',	'12.5.1',	'13.5.1',	'14.5.1',	'15.5.1',	'16.5.1',	'PN.1',	'10.5.2',	'11.5.2',	'12.5.2',	'13.5.2',	'14.5.2',	'15.5.2',	'16.5.2',	'PN.2')

library(preprocessCore)
data_matrix <- as.matrix(data_matrix)

new_matrix  <- normalize.quantiles(data_matrix, copy = TRUE)
rownames(new_matrix) <- rownames(data_matrix)
colnames(new_matrix) <- colnames(data_matrix)

write.table(new_matrix, file="", sep='\t', quote=F)
