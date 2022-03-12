plot.aki.aik.matrix <- function(AKI, AIK, bourder.out, bourder.in, colors.out, colors.in, path, foldername, subfolder, suffix, on=T){
  if (on){
    print("Plotting Aik and Aki matrice")
    pik <- ggplot2::ggplot(AIK, ggplot2::aes(target.labels, source.labels, fill=weight))+
      ggplot2::geom_raster(na.rm=T, hjust = 1, vjust = 0)+
      viridis::scale_fill_viridis(option = "A", direction = 1, na.value="white")+
      ggplot2::geom_hline(yintercept = bourder.out, color="green")+
      ggplot2::geom_vline(xintercept = bourder.out, color="green")+
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, size = 4, color=colors.out),
                     axis.text.y = ggplot2::element_text(size=4, color=colors.out %>% rev()))
    
    pki <- ggplot2::ggplot(AKI, ggplot2::aes(target.labels, source.labels, fill=weight))+
      ggplot2::geom_raster(na.rm=T, hjust = 1, vjust = 0)+
      viridis::scale_fill_viridis(option = "A", direction = 1, na.value="white")+
      ggplot2::geom_hline(yintercept = bourder.in, color="green")+
      ggplot2::geom_vline(xintercept = bourder.in, color="green")+
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, size = 4, color=colors.in),
                     axis.text.y = ggplot2::element_text(size=4, color=colors.in %>% rev()))
    
    png(sprintf("%s/%s/%s/in_matrix%s.png", path, foldername, subfolder, suffix), width = 6, height = 5, units = 'in', res = 200)
    print(pki)
    dev.off()
    
    png(sprintf("%s/%s/%s/out_matrix%s.png", path, foldername, subfolder, suffix), width = 6, height = 5, units = 'in', res = 200)
    print(pik)
    dev.off()
  } else
    print("** No Aks matrices")
}

plot.aki.aik.density<- function(AIK, AKI, path, foldername, subfolder, suffix, on=T){
  if (on){
    print("** Plotting density Aks")
    AIK$type <- "outlinks"
    AKI$type <- "inlinks"
    A <- AIK %>% rbind(AKI)
    
    ph <- ggplot2::ggplot(A, ggplot2::aes(weight, fill=type))+
      ggplot2::geom_density(alpha=0.5, na.rm = T)+
      ggplot2::xlab("similarity")+
      ggplot2::theme_classic()
    
    png(sprintf("%s/%s/%s/density%s.png", path, foldername, subfolder, suffix), width = 6, height = 5, units = 'in', res = 200)
    print(ph)
    dev.off()
  } else
    print("** No density Aks plot")
  
}

plot.aki.aik.dendrogram <- function(AKI, AIK, labels, colors.out, colors.in,  nt, path, foldername, subfolder, suffix="", on=T){
  if (on){
    print("** Plotting Aks dendrograms")
    AIK[is.na(AIK)] <- 0
    AKI[is.na(AKI)] <- 0
    
    AIK[AIK != 0] <- -log10(AIK[AIK != 0])
    AKI[AKI != 0] <- -log10(AKI[AKI != 0])
    
    AIK[AIK==0] <- max(AIK, na.rm = T) + 1
    AKI[AKI==0] <- max(AKI, na.rm = T) + 1
    
    hcluster.aik <- hclust(AIK %>% as.dist(), method = "average")
    hcluster.aik$labels <- labels
    
    
    png(sprintf("%s/%s/%s/hclust_out%s.png", path, foldername, subfolder, suffix), width = 6, height = 5, units = 'in', res = 200)
    plot(hcluster.aik %>% ape::as.phylo(), cex=0.4, tip.color = colors.out, type="fan")
    dev.off()
    
    hcluster.aki <- hclust(AKI %>% as.dist(), method = "average")
    hcluster.aki$labels <- labels[1:nt]
    
    png(sprintf("%s/%s/%s/hclust_in%s.png", path, foldername, subfolder, suffix), width = 6, height = 5, units = 'in', res = 200)
    plot(hcluster.aki %>% ape::as.phylo(), cex=0.4, tip.color = colors.in, type="fan")
    dev.off()
  } else{
    print("No Aks dendrograms")
  }
}

work.aki.aik <- function(net, nodes, leaves, mode, labels, inst, path="", foldername="", subfolder=""){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  ## Number of  columns
  nt <- 107
  suffix <- ""
  
  ## Get and format regions
  print("** Warning: Be careful selecting the right regions")
  regions <- read.csv("../CSV/Regions/Table_areas_regions_09_2019.csv")
  colnames(regions) <- c("AREA", "REGION")
  regions$AREA <- regions$AREA %>% tolower()
  
  source("functions/format_regions.R")
  regions <- format.regions(regions)
  regions <- regions[order(regions$REGION),]
  
  source("functions/get_bourder_regions.R")
  bourder.out <- get.bourder.regions(regions)
  bourder.in <- get.bourder.regions(regions[regions$AREA %in% labels[1:nt],])
  
  ## Set color for regions
  source("functions/gg_color_hue.R")
  uregions <- regions$REGION %>% unique() 
  nr <- uregions %>% length()
  color.regions <- gg.color.hue(nr)
  table.regions.out <- table(regions$REGION[regions$AREA %in% labels])
  colors.out <- c()
  for (i in 1:nr){
    colors.out <- c(colors.out, rep(color.regions[i], table.regions.out[i]))
  }
  table.regions.in <- table(regions$REGION[regions$AREA %in% labels[1:nt]])
  colors.in <- c()
  for (i in 1:nr){
    colors.in <- c(colors.in, rep(color.regions[i], table.regions.in[i]))
  }
  
  ## Get aik and aki matrices
  source("functions/compute_aki_aik.R")
  A <- compute.aki.aik(net, leaves, nodes, mode, nt, inst, suffix=suffix, save=T)
  AIK <- A$AIK
  AKI <- A$AKI
  AKI[AKI == -1] <- NA
  AIK[AIK == -1] <- NA
  
  ## Plot section 1
  colors.out.dendo <- colors.out[match(labels, regions$AREA[regions$AREA %in% labels])]
  colors.in.dendo <- colors.in[match(labels[1:nt], regions$AREA[regions$AREA %in% labels[1:nt]])]
  plot.aki.aik.dendrogram(AKI, AIK, labels, colors.out.dendo, colors.in.dendo, nt, path, foldername, subfolder, suffix, on=T)

  ## Convert AKI and AIK into data-frames with labels and ordering the labels 
  ## by regions
  source("functions/adj_to_df.R")
  AIK <- adj.to.df(AIK)
  AKI <- adj.to.df(AKI)
  AKI$source.labels <- labels[AKI$source]
  AKI$target.labels <- labels[AKI$target]
  AIK$source.labels <- labels[AIK$source]
  AIK$target.labels <- labels[AIK$target]
  
  AKI$target.labels <- factor(AKI$target.labels, levels = regions$AREA)
  AIK$target.labels <- factor(AIK$target.labels, levels = regions$AREA)
  AKI$source.labels <- factor(AKI$source.labels, levels = rev(regions$AREA))
  AIK$source.labels <- factor(AIK$source.labels, levels = rev(regions$AREA))
  
  ## Plot section 2
  plot.aki.aik.matrix(AKI, AIK, bourder.out, bourder.in, colors.out, colors.in, path, foldername, subfolder, suffix, on = T)
  plot.aki.aik.density(AIK, AKI, path, foldername, subfolder, suffix, on = T)
  
}