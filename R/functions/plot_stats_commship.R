plot.region.commship <-  function(net, hcluster, k, labels, memberships, path="", foldername="", subfolder="", filename=""){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  source("functions/assign_commship.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/format_regions.R")
  
  print("** Warning: Be careful selecting the right regions")
  regions <- read.csv("../CSV/Regions/Table_areas_regions_09_2019.csv")
  colnames(regions) <- c("AREA", "REGION")
  regions$AREA <- regions$AREA %>% tolower()
  regions <- format.regions(regions)
  
  # print("** Using hclust as passed")
  # source('functions/assign_commship.R')
  # net <- assign.commship(net, hcluster, k, 0)
  
  print("*** Warning: You are using a reference hclust to munkres the current hclust")
  # source('functions/assign_commship_reference.R')
  # net <- assign.commship.reference(net, hcluster, k)
  source('functions/format_lincomm.R')
  net <- format.lincomm(net, hcluster, k)
  
  net$commship <- net$commship %>% as.factor()
  net$source.label <- labels[net$source]
  net$target.label <- labels[net$target]

  commships <- net$commship %>% unique() %>% sort()
  nm <- commships %>% length()
  area.count <- data.frame(nodes=rep(1:length(labels), nm), regions=rep(regions$REGION[match(labels, regions$AREA)], nm), areas=rep(labels, nm), commship=rep(commships, each=length(labels)), size=0)

  for (i in 1:nrow(net)){
    area.count$size[which(area.count$nodes == net$source[i] & area.count$commship == net$commship[i])] <- area.count$size[which(area.count$nodes == net$source[i] & area.count$commship == net$commship[i])]  + 1
    area.count$size[which(area.count$nodes == net$target[i] & area.count$commship == net$commship[i])] <- area.count$size[which(area.count$nodes == net$target[i] & area.count$commship == net$commship[i])]  + 1
  }

  tr.labels <- labels[order(memberships)]
  area.count$areas <- factor(area.count$areas, levels = tr.labels)

  p <- ggplot2::ggplot(area.count, ggplot2::aes(areas, size, fill=commship))+
    ggplot2::facet_wrap(~regions, scales = "free")+
    ggplot2::geom_bar(stat = "identity")+
    ggplot2::ylab("Number of links")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90))

  png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 12, height = 5, units = 'in', res = 200)
  print(p)
  dev.off()
}

plot.wdis.commship <- function(net, hcluster, k, path="", foldername="", subfolder="", filename=""){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  source("functions/assign_commship.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  
  net <- assign.commship(net, hcluster, k, 0)
  nat <- df.to.adj(net)
  nat <- adj.to.df(nat)
  print("** Warning: Be careful choosing the right distance matrix. Be sure it has the same order as the fln matrix")
  source("functions/get_tracto2016.R")
  distances <- get.tracto2016(labels)
  distances <- distances %>% adj.to.df()
  net$commship <- net$commship %>% as.factor()
  net$distances <- distances$weight[nat$weight != 0]
  p <- ggplot2::ggplot(net, ggplot2::aes(distances, w, color=commship))+
    ggplot2::geom_point(size=0.3)+
    ggplot2::geom_smooth(method = "lm")+
    ggplot2::theme_classic()
  
  png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 6, height = 5, units = 'in', res = 200)
  print(p)
  dev.off()
}

plot.area.commship <- function(net, hcluster, k, labels, memberships, path="", foldername="", subfolder="", filename=""){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  source("functions/assign_commship.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  
  # print("** Using hclust as passed")
  # source('functions/assign_commship.R')
  # net <- assign.commship(net, hcluster, k, 0)
  
  print("*** Warning: You are using a reference hclust to munkres the current hclust")
  # source('functions/assign_commship_reference.R')
  # net <- assign.commship.reference(net, hcluster, k)
  source('functions/format_lincomm.R')
  net <- format.lincomm(net, hcluster, k)
  
  net$commship <- net$commship %>% as.factor()
  net$source.label <- labels[net$source]
  net$target.label <- labels[net$target]

  commships <- net$commship %>% unique() %>% sort()
  nm <- commships %>% length()
  area.count <- data.frame(nodes=rep(1:length(labels), nm), nodes.cluster= rep(memberships, nm), areas=rep(labels, nm), commship=rep(commships, each=length(labels)), size=0)
  
  for (i in 1:nrow(net)){
    area.count$size[which(area.count$nodes == net$source[i] & area.count$commship == net$commship[i])] <- area.count$size[which(area.count$nodes == net$source[i] & area.count$commship == net$commship[i])]  + 1
    area.count$size[which(area.count$nodes == net$target[i] & area.count$commship == net$commship[i])] <- area.count$size[which(area.count$nodes == net$target[i] & area.count$commship == net$commship[i])]  + 1
  }
  
  tr.labels <- labels[order(memberships)]
  area.count$areas <- factor(area.count$areas, levels = tr.labels)
  
  p <- ggplot2::ggplot(area.count, ggplot2::aes(areas, size, fill=commship))+
    ggplot2::facet_wrap(~nodes.cluster, scales = "free")+
    ggplot2::geom_bar(stat = "identity")+
    ggplot2::ylab("Number of links")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, size=5))
  
  png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 12, height = 5, units = 'in', res = 200)
  print(p)
  dev.off()
}