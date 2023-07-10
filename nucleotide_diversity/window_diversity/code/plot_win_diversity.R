library(scales)
library(dplyr)

# receive cmd arguments
Args <- commandArgs(T)
wd <- Args[1]
dd <- Args[2]
para <- Args[3]

# # parameters for test only
# wd <- "/Users/siyuansmac/Desktop/PhD research/Suzukii WGS analysis/manuscript_suzukii/Figures/Figure 3/plot/dot_pop_mean/plot_most_pops"
# dd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/window_diversity/results"
# para <- "winlen_125000_mincov_12_mincount_1"
# chr <- "X"

# parameter specified inside script
# selected.range <- names(vc_color)
selected.range <- c("China", "Japan", "N. America", "Europe")
# selected.range <- c("China", "Japan", "Hawaii", "N. America", "Brazil", "Europe", "Ind. O.")

# def a function to plot dots for window nucleotide diversity
plot_dot <- function(wd, dd, para, chr, selected.range) {

## read the table including Fst
setwd(wd)
df <- read.table(file=paste(dd, "/", para, "/", "diversity_", chr, "_sorted.txt", sep = ""), 
                 header=T, sep="\t", check.names=F)
df_color <- data.frame(sample=colnames(df[,-1:-3]), 
                       range=factor(c("Ind. O.", "Europe", "Europe", "China",
                                      "China", "Europe", "Europe", "Europe",
                                      "Europe", "Japan", "Japan", "China", 
                                      "N. America", "N. America", "China", "N. America", 
                                      "N. America", "N. America", "Hawaii", "N. America",
                                      "Brazil", "Europe", "N. America", "N. America",
                                      "Japan", "Europe", "Japan", "Europe", "Europe"),
                                    levels = c("Hawaii", "N. America", "Brazil", "Ind. O.", "Europe", "China", "Japan")))
df_color <- df_color[order(df_color$sample),]
df_color <- df_color[order(df_color$range),]
df_color$range <- as.character(df_color$range)

# update 10.09.2021: separate island populations from the mainland population ranges, and assign separate color to them
vc_color <- c(rgb(232, 125,	114, maxColorValue = 255), rgb(83, 182,76, maxColorValue = 255), rgb(109, 157, 248, maxColorValue = 255), rgb(109 - 35, 157 + 35, 248, maxColorValue = 255),
              rgb(232, 125,	114 + 70, maxColorValue = 255), rgb(232 -50, 125,	114 + 20, maxColorValue = 255), rgb(83 + 70, 182,76, maxColorValue = 255))
names(vc_color) <- c("N. America", "Europe", "China", "Japan", "Hawaii", "Brazil","Ind. O.")
vc_spcolor <- vc_color[df_color$range]
df_color$color <- vc_spcolor

# calculate average diversity for each geographic range
for (range in selected.range) {
  pop_to_avg <- df_color[df_color$range == range, "sample"]
  df <- df %>% mutate(!!range := rowMeans(select(., pop_to_avg)))
}

# determine coordinates of the shaded background
coord <- c(0, row.names(df[df$window == 0,][-1,]))
if (length(coord) %% 2 == 1) {
  coord <- c(coord, 1e10)
}
xleft <- coord[seq(1, length(coord), 2) + 1]
xright <- coord[seq(2, length(coord), 2) + 1]

# show the window at the middle of each contig (longer than 20 windows) that has the lowest average diversity across selected samples
mt.out <- matrix(, nrow = 0, ncol = dim(df)[2])
win.lowest.div <- c()
selected.pop <- df_color[df_color$range %in% selected.range, 1]
df.selected.pop <- select(df, one_of(selected.pop)) # select populations within the four major continental ranges
df$meandiv <- apply(df.selected.pop, 1, function(x){mean(x)})
for (contig in unique(df$contig)) {
  df.contig <- df[df$contig == contig, ]
  if (dim(df.contig)[1] >= 20) {
    df.contig.mid <- df.contig[df.contig$window > .1*length(df.contig$window) & df.contig$window < .9*length(df.contig$window), ]
    meandiv.lowest <- df.contig.mid[order(df.contig.mid$meandiv), ][1,]$meandiv
    win.lowest.div <- c(win.lowest.div,
                        rownames(df.contig.mid[df.contig.mid$meandiv <= 1.05*meandiv.lowest, ]))
    mt.out <- rbind(mt.out, df.contig.mid[df.contig.mid$meandiv <= 1.05*meandiv.lowest, ])
  }
}

# output the lowest-diversity windows to a text file
write.table(mt.out, file = paste(wd, '/', "window_lowest_diversity.txt", sep = ""),
            quote = F, sep = "\t", row.names = F, append = T, col.names = F)

## plot the dot plot
pdf(file=paste(wd, '/', "diversity_", chr, "_", para, ".pdf", sep = ""))
par(ps = 12, cex = 1, cex.main = 1)
for (range in selected.range) {
  if (range == selected.range[1]){
    plot(x = 1:nrow(df), y = df[, range], 
         col = alpha(df_color[df_color$range == range, 3], alpha = 0.7), pch = 19, cex = 0.5,
         ylim = c(0, 0.05), xlab = "Genomic window", ylab = "Nucleotide diversity", main = chr,
         panel.first = rect(xleft, 0 , xright, 0.055, col='lightgrey', border=NA),
         xaxt = "n")
    axis(1, at = round(seq(0, nrow(df), length.out = 6)))
  }
  else {
    points(x = 1:nrow(df),  y = df[, range], col = alpha(df_color[df_color$range == range, 3], alpha = 0.7), pch = 19, cex = 0.5)
  }
}

# abline(v = win.lowest.div, lty = "dotted", lwd = 0.8)
legend("topright", legend=selected.range,title="Range", cex = 0.8,
       fill=vc_color[selected.range], bty="n", xpd = T)
dev.off()
}

# plot
plot_dot(wd, dd, para, "2L", selected.range)
plot_dot(wd, dd, para, "2R", selected.range)
plot_dot(wd, dd, para, "3L", selected.range)
plot_dot(wd, dd, para, "3R", selected.range)
plot_dot(wd, dd, para, "4", selected.range)
plot_dot(wd, dd, para, "X", selected.range)
