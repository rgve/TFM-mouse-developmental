expression <- read.table('/no_backup/rg/ragarcia/repeated/MBrain/expression/excluded/selected.genes.expr.0.upreg.tsv', sep="\t" , header=F)
row.names(expression)=expression$V1
expression$V1 <- NULL

max <- apply(expression, 1, max)
min <- apply(expression, 1, min)
difference <- max-min

boxplot(difference, ylab='TPM difference', xlab='Hrt Upregulated < 1 TPM ', ylim = c(0,10) ,col = "antiquewhite",border = "chocolate")
