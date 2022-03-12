plot.aks.dendrogram <- function(AKS, labels, regions, colors, path="", foldername="", subfolder="", preffix="", suffix=""){
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  colors.dendo <- colors[match(labels, regions$AREA[regions$AREA %in% labels])]
  
  AKS[is.na(AKS)] <- 0
  AKS[AKS != 0] <- -log10(AKS[AKS != 0])
  AKS[AKS==0] <- max(AKS, na.rm = T) + 1
  
  hcluster.aks <- hclust(AKS %>% as.dist(), method = "average")
  hcluster.aks$labels <- labels
  
  png(sprintf("%s/%s/%s/%s_hclust%s.png", path, foldername, subfolder, preffix, suffix), width = 6, height = 5, units = 'in', res = 200)
  plot(hcluster.aks %>% ape::as.phylo(), cex=0.4, tip.color = colors.dendo, type="fan")
    dev.off()
}