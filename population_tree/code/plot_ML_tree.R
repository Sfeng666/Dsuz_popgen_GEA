# 1. receive cmd arguments
Args <- commandArgs(T)
wd <- Args[1]
dd <- Args[2]
sp <- Args[3]
cl <- Args[4]
para <- Args[5]

# parameters for test only
wd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/ML_tree/plot/mincount_10_no_se_auto_X_over10"
dd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/ML_tree/results/mincount_10_no_se_auto_X_over10"
sp <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/ML_tree/treemix/treemix-1.13/src/plotting_funcs.R"
cl <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/ML_tree/plot/pop_order_color_consistent.txt"
para <- "winlen_500_mincov_12_mincount_10"
# chr <- "X"

# 2. load libraries
source(sp)
library(zoo)

# 3. def a function to plot ML tree, residual tree and likelihood line
plot.MLtree <- function(wd, dd, sp, cl, para, chr) {
  
## 3.1. create a sub directory
dir.create(file.path(wd, para), showWarnings = FALSE)
setwd(file.path(wd, para))

## 3.2. plot ML tree and residual matrix for each number of migration events
for (m in 0:20){
       treeout <- paste(dd, "/", para, "/", chr, "_m_", m, sep = "")
       pdf(file = paste(wd, "/", para, "/", chr, "_m_", m, '.pdf', sep = ""), width = 20, height = 10)
       par(mfrow = c(1,2), ps = 11, cex = 1, cex.main = 1)
       plot_tree(treeout, cl, disp = 0.00002, plus = 0.002)
       plot_resid(treeout, cl)
       title(main = paste(para, '_mig_', m, sep = ""), line = -3, outer = TRUE)
       dev.off()
}

## 3.3. plot a curve of likelihood based on the number of migration events
### 3.3.1. get the ln(likelihood) transformed likelihood value
ln.llik <- sapply(0:20, function(m){
  fileName <- paste(dd, "/", para, "/", chr, "_m_", m, ".llik", sep = "")
  line.llik <- readLines(file(fileName,open="r"))[2]
  llik <- as.numeric(gsub("^.*: (-[0-9]+.*$)", "\\1", line.llik))
})

### 3.3.2. find approximately where the curve reaches plateau
x <- 0:20
y <- ln.llik

slope <- diff(y) # Calculate the slope of the curve
smooth_slope <- rollmean(slope, k=5, align="right") # Smooth the slope data using a moving average

#### Identify the point where the slope becomes approximately zero
threshold <- 0.1 * max(abs(smooth_slope)) # set the plateu threshold as 10% of the maximum smoothed slope
plateau <- min(which(abs(smooth_slope) < threshold))
x_plateau <- x[plateau] # Find the corresponding x-value for the identified plateau

### 3.3.3. plot the curve
pdf(file = paste(wd, "/", para, "/", chr, "_likelihood", '.pdf', sep = ""))
par(ps = 11, cex = 1, cex.main = 1)
plot(0:20, ln.llik, type="l", lwd=2, 
     xlab="Number of migration events", ylab="ln(likelihood)", main=paste(chr, "-linked contigs", sep = ""),
     yaxt="n")
axis(2, at=pretty(ln.llik, n=5)) # create y axis with 5 break points
abline(v=x[plateau], col="red")
dev.off()

## 3.4. plot a curve of the fraction of explained variance based on the number of migration events
### 3.4.1. define a function to calculate f
calc_f <- function(dir.cov, file.cov.est, file.cov.fitted){
  path.cov.est <- gzfile(paste(dir.cov, file.cov.est, sep = '/'))
  cov.est <- as.matrix(read.table(path.cov.est, header=T))
  cov.est <- cov.est[order(rownames(cov.est)), order(colnames(cov.est))] # the fitted covariance matrix is in a different order than the sample covariance matrix. Symmetrical reordering needed!
  
  path.cov.fitted <- gzfile(paste(dir.cov, file.cov.fitted, sep = '/'))
  cov.fitted <- as.matrix(read.table(path.cov.fitted, header=T))
  cov.fitted <- cov.fitted[order(rownames(cov.fitted)), order(colnames(cov.fitted))] # the fitted covariance matrix is in a different order than the sample covariance matrix. Symmetrical reordering needed!
  
  cov.est.upper <- cov.est[upper.tri(cov.est)]
  cov.est.upper.mean <- mean(cov.est.upper)
  R.upper <- (cov.est - cov.fitted)[upper.tri((cov.est - cov.fitted))]
  R.upper.mean <- mean(R.upper)
  f <- 1 - sum((R.upper - R.upper.mean)^2)/sum((cov.est.upper - cov.est.upper.mean)^2)
  return(f)
}

### 3.4.2. get the value of  the fraction of explained variance
v.frac <- sapply(0:20, function(m){
  dir.cov <- paste(dd, "/", para, sep = "")
  file.cov.est <- paste(chr, "_m_", m, ".cov.gz", sep = "")
  file.cov.fitted <- paste(chr, "_m_", m, ".modelcov.gz", sep = "")
  calc_f(dir.cov, file.cov.est, file.cov.fitted)
})

### 3.3.3. find approximately where the curve reaches plateau
x <- 0:20
y <- v.frac

slope <- diff(y) # Calculate the slope of the curve
smooth_slope <- rollmean(slope, k=5, align="right") # Smooth the slope data using a moving average

#### Identify the point where the slope becomes approximately zero
threshold <- 0.1 * max(abs(smooth_slope)) # set the plateau threshold as 10% of the maximum smoothed slope
plateau <- min(which(abs(smooth_slope) < threshold))
x_plateau <- x[plateau] # Find the corresponding x-value for the identified plateau

### 3.4.4. plot the curve
pdf(file = paste(wd, "/", para, "/", chr, "_fraction_explained_var", '.pdf', sep = ""))
par(ps = 11, cex = 1, cex.main = 1)
plot(0:20, v.frac, type="l", lwd=2, 
     xlab="Number of migration events", ylab="Fraction of variance explained by the model", main=paste(chr, "-linked contigs", sep = ""),
     yaxt="n", ylim = c(min(v.frac), 1))
axis(2, at=pretty(v.frac, n=4)) # create y axis with 5 break points
abline(v=x[plateau], col="red")
dev.off()
}

# 4. plot
plot.MLtree(wd, dd, sp, cl, para, "auto")
plot.MLtree(wd, dd, sp, cl, para, "X")