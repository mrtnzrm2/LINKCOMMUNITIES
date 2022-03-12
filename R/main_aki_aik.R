make.aks.matrix <- function(AKS, bourder, colors, labels, regions, inst, preffix="", suffix="", on=T){
  if (on){
    print("** Plotting Aks matrix")
    source("functions/plot_aks_matrix.R")
    plot.aks.matrix(AKS, bourder, colors, labels, regions, path=inst$plot, foldername=inst$mainfolder, subfolder="Preprocessing", preffix=preffix, suffix=suffix)
  } else
    print("** No Aks matrix")
}

make.aks.dendrogram <- function(AKS, labels, regions, colors, inst, preffix="", suffix="", on=T){
  if (on){
    print("** Plotting Aks dendrogram")
    source("functions/plot_aks_dendrogram.R")
    plot.aks.dendrogram(AKS, labels, regions, colors, path=inst$plot, foldername=inst$mainfolder, subfolder="Preprocessing", preffix=preffix, suffix=suffix)
  } else
    print("** No Aks dendrogram")
}

make.aki.aik.density <- function(AKI, AIK, labels, regions, inst, suffix="", on=T){
  if (on){
    print("** Plotting Aks density")
    source("functions/plot_aki_aik_density.R")
    plot.aki.aik.density(AKI, AIK, labels, regions, path=inst$plot, foldername=inst$mainfolder, subfolder="Preprocessing", suffix=suffix)
  } else
    print("** No Aks density")
}

make.aki.aik.scatter <- function(net, AKI, AIK, nt, labels, inst, suffix="", on=T){
  if (on){
    print("** Plotting scatter plots")
    source("functions/plot_aki_aik_scatter.R")
    plot.aki.aik.scatter(net, AKI, AIK, nt, labels, path=inst$plot, foldername=inst$mainfolder, subfolder="Preprocessing", suffix=suffix)
  } else
    print("** No scatter plots")
}

main <- function(inst){
  library(magrittr)
  
  linkage <- 'average'
  index <- 'jaccp'
  mode <- 'ALPHA'
  nlog10 <- T
  tag <- paste(inst$model, inst$distances, inst$folder, sep = "_")
  fln.name <- 'fln'
  nt <- 107   ## Number of  columns
  suffix <- ""
  
  source('functions/model_name_analysis.R')
  filename <- model.name(index, linkage, mode, nlog10, F, tag)
  source('functions/sformat.R')
  inst$mainfolder <- paste(inst$folder, inst$common, paste(toupper(linkage), 'full','l10', sep = "_"), sep = "/")
  dir.create(sprintf('%s/%s', inst$plot, inst$mainfolder), showWarnings = F)
  
  #### Load network
  source('functions/load_net.R')
  netx <- load.net(inst)
  net <- netx$net
  nodes <- netx$nodes
  leaves <- netx$leaves
  labels <- netx$labels
  source("functions/format_labels.R")
  labels <- labels %>% format.labels()
  coords <- netx$fln.coords
  regions <- netx$regions
  
  if (nlog10){
    net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
  }
  
  source("functions/get_regions_colors.R")
  cr <- get.regions.colors(labels, nt)
  
  ## Get aik and aki matrices
  source("functions/compute_aki_aik.R")
  A <- compute.aki.aik(net, nodes, mode, nt, inst, suffix=suffix, save=T)
  AIK <- A$AIK
  AKI <- A$AKI
  AKI[AKI == -1] <- NA
  AIK[AIK == -1] <- NA
  
  ## Plotting
  make.aks.dendrogram(AKI, labels[1:nt], cr$rgn, cr$cin, inst, preffix="in", suffix=suffix, on=T)
  make.aks.dendrogram(AIK, labels,cr$rgn, cr$cout, inst, preffix="out", suffix=suffix, on=T)
  make.aks.matrix(AKI, cr$bin, cr$cin, labels, cr$rgn, inst, preffix="in", suffix=suffix, on = T)
  make.aks.matrix(AIK, cr$bout, cr$cout, labels, cr$rgn, inst, preffix="out", suffix=suffix, on = T)
  make.aki.aik.density(AKI, AIK, labels, cr$rgn, inst, suffix=suffix, on=T)
  make.aki.aik.scatter(net, AKI, AIK, nt, labels, inst, suffix=suffix, on=T)

}


### HEAD ###
folder <- 'merged'
distances <- 'tracto2016'
model <- 'zz_model'
csv.name <- 'fln'
labels.name <- 'imputation_labels'  # imputation_labels  rlabels

common.path <- paste(distances, model, sep = '/')
csv.path <- paste(folder, 'imputation', common.path, paste0(csv.name,'.csv'), sep = '/') # merged/imputation/tracto2016/zz_model/   91x40
labels.path <- paste(folder, 'labels', common.path, paste0(labels.name,'.csv'), sep = '/')
plot.path <- '../plots'

path.list <- list(csv=csv.path, 
                  labels=labels.path, 
                  plot=plot.path, 
                  common=common.path,
                  folder=folder,
                  model=model,
                  distances=distances)

main(path.list)