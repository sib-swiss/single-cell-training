
plotMeanVar2 = function( count_mat ){
  log_mean <- log10(rowMeans(count_mat))
  log_var <- log10(apply(count_mat, 1, var))
  cols <- densCols(log_mean, log_var)
  plot(log_mean, log_var, col = cols, pch = 19, cex = 0.2, xlab = "log10(mean)", ylab = "log10(variation)")
  abline(0,1, col = "red")
}

for(type in c("counts", "data", "scale.data")){
  png(paste0("plot_mean_var_", type, ".png"), pointsize = 14)
  plotMeanVar2(GetAssayData(gbm, slot =  type))
  dev.off()
}

