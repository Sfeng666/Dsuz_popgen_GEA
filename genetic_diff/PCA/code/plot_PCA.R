library(stats)
library(data.table)
library(readr)
library("scatterplot3d")
library(robust)

# receive cmd arguments
Args <- commandArgs(T)
wd <- Args[1]
dd <- Args[2]
para <- Args[3]

# # parameters for test only
# wd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/PCA/plot"
# dd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/PCA/results"
# para <- "mincov_12_mincount_1"
# contig <- "auto"

scatter_PCA <- function(wd, dd, para, contig) {
  
setwd(wd)
# df <- read.table(file=paste(dd, "/", para, "/", "genetic_distance_", contig, ".txt", sep = ""), 
#                  header=T, sep="\t", check.names=FALSE)
df <- fread(file=paste(dd, "/", para, "/", "maf_", contig, ".txt.gz", sep = ""),
            header = T, showProgress = T, data.table = F, sep="\t", check.names=FALSE)

# calculate PCs
pc <- prcomp(t(df),center=TRUE)

# summarize the proportion of the top three PC
pc_sum <- summary(pc)
proportion_pc1 <- pc_sum$importance[2,1]
proportion_pc2 <- pc_sum$importance[2,2]
proportion_pc3 <- pc_sum$importance[2,3]

# add color to population ranges
pc_x <- as.data.frame(pc_sum$x)
pc_x$range <- c("FR-Run", "European", "European", "Chinese",
               "Chinese", "European", "European", "European",
               "European", "Japanese", "Japanese", "Chinese", 
               "American", "American", "Chinese", "American", 
               "American", "American", "US-Haw", "American",
               "BR-Pal", "European", "American", "American",
               "Japanese", "European", "Japanese", "European", "European")
vc_color <- c(rgb(232, 125,	114, maxColorValue = 255), rgb(83, 182,76, maxColorValue = 255), rgb(109, 157, 248, maxColorValue = 255), rgb(109 - 35, 157 + 35, 248, maxColorValue = 255),
              rgb(232, 125,	114 + 70, maxColorValue = 255), rgb(232 -50, 125,	114 + 20, maxColorValue = 255), rgb(83 + 70, 182,76, maxColorValue = 255))
vc_color_border <- c(rgb(232/2, 125/2, 114/2, maxColorValue = 255), 
                     rgb(83/2, 182/2, 76/2,	maxColorValue = 255), 
                     rgb(109/2, 157/2, 248/2, maxColorValue = 255),
                     rgb((109 - 35)/2, (157 + 35)/2, 248, maxColorValue = 255),
                     rgb(232/2, 125/2, (114 + 70)/2, maxColorValue = 255),
                     rgb((232 - 50)/2, 125/2, (114 + 20)/2, maxColorValue = 255),
                     rgb((83 + 70)/2, 182/2, 76/2, maxColorValue = 255))
names(vc_color) <- c("American", "European", "Chinese", "Japanese", "US-Haw", "BR-Pal","FR-Run")
names(vc_color_border) <- c("American", "European", "Chinese", "Japanese", "US-Haw", "BR-Pal","FR-Run")
vc_spcolor <- vc_color[pc_x$range]
vc_spcolor_border <- vc_color_border[pc_x$range]
pc_x$color <- vc_spcolor
pc_x$color_border <- vc_spcolor_border

# identify outliers by calculating Robust Mahalanobis distance on the top three PCs
# Note: Mahalanobis distance is a multivariate distance based on all variables (PCs here) at once. 
# Ref: https://www.r-bloggers.com/2019/08/detecting-outlier-samples-in-pca/
dist <- covRob(pc_x[,1:3], estim = "pairwiseGK")$dist
outliers <- names(dist[dist > quantile(dist, probs = 0.90)])

## plot PCA
pdf(file=paste(wd, '/', "PCA_", contig, "_", para, ".pdf", sep = ""))
par(ps = 12, cex = 1, cex.main = 1)
s3d <- scatterplot3d(x = pc_x$PC1, y = pc_x$PC2, z = pc_x$PC3, 
                     type = "h",
                     color = pc_x$color_border,
                     cex.symbols = 2,
                     xlab = paste('PC1 (explained variance, ', round(proportion_pc1*100, 2), '%)', sep = ""), 
                     ylab = paste('PC2 (explained variance, ', round(proportion_pc2*100, 2), '%)', sep = ""),
                     zlab = paste('PC3 (explained variance, ', round(proportion_pc3*100, 2), '%)', sep = ""),
                     bg = pc_x$color,
                     pch = 21, 
                     grid = T,
                     main = para)
legend("topleft", title = 'Range', names(vc_color), fill=vc_color, bty="n")

# label outliers
text(s3d$xyz.convert(pc_x[rownames(pc_x) %in% outliers, 1:3]), 
     labels = outliers,
     cex= 0.7, 
     adj = c(-0.1, -1),
     col = rgb(37,38,40, maxColorValue = 255))

dev.off()

# remove the variable and clear memory
rm(df)
gc()
}

# plot
scatter_PCA(wd, dd, para, "auto")
scatter_PCA(wd, dd, para, "X")