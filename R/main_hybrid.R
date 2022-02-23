load.hclust <- function(resolution, filename, suffix, save=T){
  
  if (save){
    merge <- read.csv(paste('..', 'CSV', 'merged/similarity/tracto2016/zz_model', paste0('merge_th_', resolution, paste0(suffix,'.csv')), sep = '/'), header = F) %>% as.matrix()  # 91x40
    height <- read.csv(paste('..', 'CSV', 'merged/similarity/tracto2016/zz_modeL', paste0('height_th_', resolution, paste0(suffix,'.csv')), sep = '/'), header = F) %>% as.matrix()
    
    hybrid <- list(merge = merge,
                   height = height,
                   order = 1:leaves)
    class(hybrid) <- 'hclust'
    
    saveRDS(hybrid, paste('../RDS', filename,'hclust', paste0('hierarchical_clustering_hyb_', resolution,'.rds'), sep = '/')) ## save
  } else{
    hybrid <- readRDS(paste('../RDS', filename,'hclust', paste0('hierarchical_clustering_hyb_', resolution,'.rds'), sep = '/')) 
  }
  
  return(hybrid)
}

load.process <- function(net, hybrid, nodes, filename, save=T){
  if (save){
    source('functions/process_hclust.R')
    hybrid.process <- process.hclust(net, hybrid, nodes)
    saveRDS(hybrid.process, file.path('..', 'RDS', '{name}/dndrgm/single_linkage_features_hyb.rds' %>% sformat(list(name=filename))) )
  }else{
    hclust.features <- readRDS(file.path('..', 'RDS', '{name}/dndrgm/single_linkage_features_hyb.rds' %>% sformat(list(name=filename))) )
  }
  
  return(hclust.features)
}

load.merde <- function(process, nodes, labels, net, hybrid, leaves, resolution, filename, save=T){
  if (save){
    source('functions/la_arbre_a_merde.R')
    hybrid.merde <- la.arbre.a.merde(process, nodes, labels[1:nodes], net, hybrid, leaves)
    saveRDS(hybrid.merde, paste('../RDS', filename, 'merde', paste0('node_agglomeration_hyb_', resolution,'.rds'), sep = '/')) 
  } else{
    hybrid.merde <- readRDS(paste('../RDS', filename, 'merde', paste0('node_agglomeration_hyb_', resolution,'.rds'), sep = '/'))
  }
  return(hybrid.merde)
}

make.parameters <- function(inst, folder.plot, features, on=T){
  if (on){
    print("Printing parameters plots")
    source('functions/plot_process_parameters.R')
    plot.process.parameters(inst$plot, folder.plot, features)  
  } else
    print("No parameters")
}

make.link.comm <- function(Ns, net, nodes, cluster, leaves, regions, inst, folder.plot, on=T){
  if (on){
    print("Plotting link communities matrices")
    source('functions/plot_lincomm_matrix.R')
    for (k.2 in Ns){
      plot.lincomm.matrix(net[net$source <= nodes & net$target <= nodes,], cluster,
                          leaves, nodes, k.2, 0,
                          labels, regions, plt=T, subfolder='Matrices',
                          path=inst$plot,
                          filename = sprintf('%s', k.2),
                          foldername = folder.plot)
    } 
  } else
    print("No link communities matrix")
}

save.tables <- function(Ks, net, merde, labels, inst, resolution, suffix, on=T){
  if (on){
    print("Creating tables")
    source("functions/link_table.R")
    for (k in Ks){
      memberships <- link.table(net, merde, labels, k)
      h.memberships <- memberships[memberships$DIC == 'H',]
      h.memberships <- h.memberships[match(labels, h.memberships$AREA),]
      v.memberships <- memberships[memberships$DIC == 'V',]
      v.memberships <- v.memberships[match(labels, v.memberships$AREA),]
      dir.create(paste('../CSV', inst$folder, 'tables', inst$common, paste0('TH_', resolution, suffix), k, sep = '/'), showWarnings = F)
      write.csv(h.memberships, paste('../CSV', inst$folder, 'tables', inst$common, paste0('TH_', resolution, suffix), k,'h.csv', sep = '/'), row.names = F)
      write.csv(v.memberships, paste('../CSV', inst$folder, 'tables', inst$common, paste0('TH_', resolution, suffix), k,'v.csv', sep = '/'), row.names = F)
    } 
  } else
    print("No creating tables")
}


main <- function(inst){
  
  linkage <- 'single'
  similarity.index <- 'jaccp'
  self.loop.type <- 'ALPHA'
  EC <- F
  nlog10 <- T
  tag <- 'zz_model_NONULL_tracto2016_merged' # zz_model_NONULL_tracto2016_merged  91x40
  fln.name <- 'fln'
  single.resolution <- 2351
  suffix <- '_l10'
  
  source('functions/load_net.R')
  netx <- load.net(inst)
  net <- netx$net
  leaves <- netx$leaves
  nodes <- netx$nodes
  labels <- netx$labels
  regions <- netx$regions
  
  if (nlog10){
    net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
  }
  
  source('functions/model_name_analysis.R')
  filename <- model.name(similarity.index, linkage, self.loop.type,  nlog10, EC, tag)
  source('functions/sformat.R')
  folder.plot <- paste(inst$plot, inst$folder, inst$common, paste0('normal_full',suffix), sep = "/")
  dir.create(sprintf('%s/%s', inst$plot, folder.plot), showWarnings = F)
  
  #### Analysis stars here:
  hybrid <- load.hclust(single.resolution, filename, suffix, save=F)
  process <- load.process(net, hybrid, nodes, filename, save=F)
  print(process)
  hybrid.merde <- load.merde(process, nodes, labels, net, hybrid, leaves, single.resolution, filename, save=F)
  
  ### Plotting part:
  make.parameters(inst, folder.plot, process, on=F)
  make.link.comm(30:2, net, nodes, hybrid, leaves, regions, inst, folder.plot, on=F)
  ## Table part:
  save.tables(Ks, net, hybrid.merde, labels, inst, single.resolution, suffix, on=F)
  
}

folder <- 'merged'
distances <- 'tracto2016'
model <- 'zz_model'
csv.name <- 'fln'
labels.name <- 'imputation_labels'

common.path <- paste(distances, model, sep = '/')
csv.path <- paste(folder, 'imputation', common.path, paste0(csv.name,'.csv'), sep = '/') # merged/imputation/tracto2016/zz_model/   91x40
labels.path <- paste(folder, 'labels', common.path, paste0(labels.name,'.csv'), sep = '/')
plot.path <- '../plots'

path.list <- list(csv=csv.path, 
                  labels=labels.path, 
                  plot=plot.path, 
                  common=common.path,
                  folder=folder)

main(path.list)

# ### FOR IMPUTATION and using partional matrix
# get.edc <- net$source <= nodes & net$target <= nodes
# leaves <- net[get.edc, ] %>% nrow()
# net$id[get.edc] <- 1:leaves
# net$id[!get.edc] <- (leaves+1):nrow(net)
# net <- net[order(net$id),]
# ###

# net.sim <- readRDS(sprintf('../RDS/%s/similarity/link_similarity.rds', filename))
# net.dis <- -log10(net.sim)
# net.dis[is.na(net.dis)] <- max(net.dis, na.rm = T) + 1

