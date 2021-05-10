library(ggplot2)
library(cowplot)
getwd()
setwd("/no_backup/rg/ragarcia/MBrain")

function1 <- function(m) {
  m <- as.data.frame(t(m))
  m.scaled <- as.data.frame(scale(m))
  m.scaled <- as.data.frame(t(na.omit(t(m.scaled))))
  return(m.scaled)
} 

# read expression matrix
x <- list()
x[[1]] <- read.table("expression/MBrain/selected.genes.rep.1.2.log.after.QN.merged2.tsv", h=T, sep
                     ="\t")

# the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H3K36me3", "H3K4me3", "H3K4me1", "H3K4me2", "H3K9me3", "H3K27me3")

# read marks dfs
for ( i in 1:8 ) {
  x[[(i+1)]] <- read.table(paste0(marks[i], "/", marks[i], ".matrix.after.QN.merged.tsv"),
                           h=T, sep="\t")
}


#removing 10.5 tp
#x <- lapply(x, function(x) x[!names(x) %in% c('C10.1', 'X10.5')])

# check rownames order
x[[1]] <- x[[1]][rownames(x[[2]]), ]
for ( i in 2:9 ) {
  
  stopifnot(identical(rownames(x[[1]]),
                      rownames(x[[i]])))
}


# 6. get transposed & scaled matricesfor each mark

x.scaled <- list()
for ( i in 1:9 ) {
  
  x.scaled[[i]] <- function1(m = x[[i]])
  
}

# 7. obtain genes present in the 10 matrices
my.genes <- Reduce(intersect, list(colnames(x.scaled[[1]]),
                                   colnames(x.scaled[[2]]),
                                   colnames(x.scaled[[3]]),
                                   colnames(x.scaled[[4]]),
                                   colnames(x.scaled[[5]]),
                                   colnames(x.scaled[[6]]),
                                   colnames(x.scaled[[7]]),
                                   colnames(x.scaled[[8]]),
                                   colnames(x.scaled[[9]])))
df <- data.frame(stringsAsFactors = F)
for ( i in 1:9 ) {
  
  df <- rbind(df, x.scaled[[i]][, my.genes])
}
pca.res <- prcomp(df, center = F, scale = F)
summary(pca.res)

# 9. prepare df for plot
y <- as.data.frame(pca.res$x[, 1:4])
y$tp_type <- c(paste0(rownames(y)[1:8], "_expression"),
               paste0(rownames(y)[1:8], "_H3K27ac"),
               paste0(rownames(y)[1:7], "_H3K9ac"),
               paste0(rownames(y)[1:8], "_H3K4me3"),
               paste0(rownames(y)[1:8], "_H3K4me1"),
               paste0(rownames(y)[1:8], "_H3K36me3"),
               paste0(rownames(y)[1:7], "_H3K4me2"),
               paste0(rownames(y)[1:8], "_H3K9me3"),
               paste0(rownames(y)[1:8], "_H3K27me3"))
y$tp <- rep(rownames(y))
y$type <- rep(c("expression", marks), times = c(8,8,7,8,8,8,7,8,8))
y$type <- factor(y$type, levels = c("expression", marks))
#y$tp_numeric <- rep(c(0,3,6,9,12,18,24,36,48,72,120,168), 9)
y$tp_numeric <- c('X10.5','X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN','X10.5','X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN', 'X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN','X10.5', 'X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN','X10.5','X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN', 'X10.5','X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN', 'X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN','X10.5', 'X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN','X10.5','X11.5', 'X12.5', 'X13.5', 'X14.5', 'X15.5', 'X16.5', 'XPN')
#y$tp_cat <- c(1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7)
y$tp_cat <- c(1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8)


# 10. define color palette
palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00",
             "expression" = "white")

palette2 <- c("H3K9ac" = "#af4d85",
              "H3K27ac" = "#630039",
              "H3K4me3" = "#d199b9",
              "H3K27me3" = "#1d2976",
              "H3K9me3" = "#a7add4",
              "H3K36me3" = "#7fbc41",
              "H3K4me1" = "#e5ab00",
              "H3K4me2" = "#a67c00",
              "expression" = "gray")

# 11. make plot
p <- ggplot(y,
            aes(x=PC1, y=PC2, fill = type, label = tp_numeric, colour=type)) +
  theme_bw() +
  facet_wrap(~type, nrow=2) +
  theme(axis.title = element_text(size =18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 90, vjust = .5, hjust = .95),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background.x = element_blank(),
        strip.text.x = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.position = "bottom",
        plot.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  geom_point(data = y,
             aes(size=tp_numeric), shape=21, alpha = .8, colour="black") +
  geom_line(data = y, alpha=.3, size=10) +
  scale_fill_manual(values = palette) +
  scale_colour_manual(values = palette2) +
  #  xlab("PC1 (13.8 %)") +
  #  ylab("PC2 (5.1 %)") +
  xlab(summary(pca.res)$importance[2,1]*100) +
  ylab(summary(pca.res)$importance[2,2]*100) +
  guides(fill=F,
         colour=F,
         size = guide_legend(title = "hours", nrow = 1)) 
# + scale_size_continuous(breaks = c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168))

my.legend <- get_legend(p)
p <- p + guides(size=F)

plot_grid(plotlist = list(p, my.legend),
          nrow = 2, rel_heights = c(1, 0.1))
