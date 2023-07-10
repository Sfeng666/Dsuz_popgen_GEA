library(ggplot2)
library(dplyr)

# def a function to format decimals in ggplot
fmt_dcimals <- function(decimals=2){
  function(x) format(x, nsmall = decimals,scientific = FALSE)
}

# read the table including divergence
setwd("/Users/siyuansmac/Desktop/PhD research/Suzukii WGS analysis/dataset.nosync/annotation/Refseq/reannotate/calc_divergence/plot")
df <- data.frame(t(read.delim(file="../divergence.txt", header=F, sep="\t", check.names=FALSE)))
df <- rename(df, cate = X1, div = X2)
df[,2] <- data.frame(sapply(df[,2], function(x) as.numeric(as.character(x))))
df$cate <- factor(df$cate, levels = df$cate)

# barplot
ggplot(data = df, mapping = aes(x = cate, y = div)) + geom_bar(stat = 'identity') +
  xlab("Category of reannotations") +
  ylab("Divergence Suz vs Bia") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.ticks = element_line(colour = "black"))
ggsave('divergence_categories.pdf')