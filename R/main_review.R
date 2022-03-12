load.merde <- function(process, nodes, labels, net, net.cluster, leaves, filename, save=T){
  if (save){
    print("Loading the merde")
    dir.create(sprintf('../RDS/%s/merde/', filename), showWarnings = F)
    source('functions/la_arbre_a_merde.R')
    ham.merde <- la.arbre.a.merde(process, nodes, labels[1:nodes], net, net.cluster, leaves)
    saveRDS(ham.merde,  paste('../RDS', filename, 'merde', 'node_agglomeration.rds', sep = '/'))
  } else{
    print("Not loading the merde")
    ham.merde <- readRDS(paste('../RDS', filename, 'merde', 'node_agglomeration.rds', sep = '/'))
  }
  return(ham.merde)
}

make.cascade <- function(net, merde, nodes, labels, regions, inst, on=T){
  if (on){
    print("Plotting the merde cascade")
    source('functions/cascade_sil_vous_plaît.R')
    cascade.merde <- cascade.sil.vous.plaît(net, merde, nodes)
    source('functions/cascade_arc_en_ciel.R')
    cascade.arc.en.ciel(merde, nodes, labels[1:nodes], regions, cascade.merde, path=inst$plot,
                        subfolder='ArcInCiel', plt=T, filename = 'Arco_byH',
                        foldername = inst$mainfolder) 
  } else
    print("No cascade")
}

make.parameters <- function(inst, features, on=T){
  if (on){
    print("Printing parameters plots")
    source('functions/plot_process_parameters.R')
    plot.process.parameters(inst$plot, inst$mainfolder, features)  
  } else
    print("No parameters")
}

make.mountain <- function(merde, labels, nodes, coords, regions, inst, on=T){
  if (on){
    print("Plotting the mountain view")
    source('functions/mountain_view.R')
    mountain.view(merde,
                  labels = labels[1:nodes],
                  net.coord = coords,
                  regions = regions,
                  path = inst$plot,
                  nodes = nodes,
                  subfolder = 'Mountain',
                  filename='ContourPlot',
                  foldername = inst$mainfolder,
                  plt=T) 
  } else
    print("No mountain view")
  
}

make.circular <- function(Ns, features, merde, inst, regions, on=T){
  if (on){
    print("Plotting circular dendrogram")
    source('functions/plot_dendrogram_circular.R')
    for (k.m.2 in Ns){
      k.2 <- features$K[which(features$NEC == k.m.2)]
      if (length(k.2) == 0)
        next
      k.2 <- k.2[length(k.2)]
      plot.dendrogram.circular(merde, k=k.m.2,
                               filename = sprintf('polarMerde_k_%i', k.m.2),
                               plt = T, alternative = F,
                               tip.labels.size =3.5,
                               line.size=1,
                               legend.plot.size = 30, strip.text.size = 1,
                               offset.strip.text = 0.05, offset.strip = 0.5,
                               branch.none = T,
                               polar.text.size= 5, bar.line.size = 3,
                               subfolder = 'BrainEvolutionSim',
                               path=inst$plot,
                               tree.regions = list(T, regions),
                               foldername = inst$mainfolder)
    }
  } else
    print("No circular dendrogram")
}

make.dend.matrix <- function(Ns, features, merde, cluster, nodes, inst, regions, on=T){
  if (on){
    print("Plotting matrix dendrogram")
    source('functions/plot_dendrogram_matrix.R')
    for (k.m.2 in Ns){
      k.2 <- features$K[which(features$NAC == k.m.2)]
      if (length(k.2) == 0)
        next
      k.2 <- k.2[length(k.2)]
      plot.dendrogram.matrix(net[net$source <= nodes & net$target <= nodes,], merde, cluster, k=k.2, height = 0,
                             labels[1:nodes], regions, nodes, k.merde=k.m.2,
                             filename = sprintf('recMerde_k_%i_v2', k.m.2),
                             foldername = inst$mainfolder,
                             plt = T, alternative = F,
                             line.size=1,
                             tip.labels.size=1,
                             branch.none = F,
                             subfolder = 'Matrices_ns',
                             path=inst$plot)
    }
  } else
    print("No matrix dendrogram")
}

make.link.comm <- function(Ns, net, nodes, labels, cluster, leaves, regions, memberships, inst, on=T){
  if (on){
    print("Plotting link communities matrices")
    source('functions/plot_lincomm_matrix.R')
    for (k.2 in Ns){
      plot.lincomm.matrix(net[net$source <= nodes & net$target <= nodes,], cluster,
                          leaves, nodes, k.2, 0,
                          labels, regions, memberships=memberships,
                          plt=T, subfolder='Matrices',
                          path=inst$plot,
                          filename = sprintf('%s_linkcommunities_%s', k.2, inst$suffix),
                          foldername = inst$mainfolder)
    } 
  } else
    print("No link communities matrix")
}

make.bars <-  function(Ks, net, merde, cluster, labels, regions, inst, on=T){
  if (on){
    print("Plotting link comms bar plots")
    source("functions/plot_bar_link_comms.R")
    for (k in Ks)
      plot.bar.link.comms(k, net, merde, cluster, labels, regions, path = inst$plot, foldername = inst$mainfolder, subfolder = "BarLinkComms", filename = paste0(k,"_bar"))  
  } else
    print("No bars")
  
}

save.tables <- function(Ks, net, merde, labels, inst, on=T){
  if (on){
    print("Creating tables")
    source("functions/link_table.R")
    for (k in Ks){
      dir.create(paste('../CSV', inst$folder, 'tables', inst$common, k, sep = '/'), showWarnings = F)
      memberships <- link.table(net, merde, labels, k)
      h1 <- memberships[memberships$DIC == 'H',]
      v1 <- memberships[memberships$DIC == 'V',]
      commships <- memberships$COMMSHIP %>% unique() %>% sort()
      nc <- commships %>% length()
      rc <- 1:nc
      for (j in 1:nc){
        h2 <- h1[h1$COMMSHIP == commships[j],]
        if (nrow(h2) < length(labels)){
          h2 <- h2 %>% rbind(data.frame(AREA=labels[!labels %in% h2$AREA], COMMSHIP=commships[j], P=0, DIC="H"))
        }
        h2 <- h2[match(labels, h2$AREA),]
        v2 <- v1[v1$COMMSHIP == commships[j],]
        if (nrow(v2) < length(labels)){
          v2 <- v2 %>% rbind(data.frame(AREA=labels[!labels %in% v2$AREA], COMMSHIP=commships[j], P=0, DIC="V"))
        }
        v2 <- v2[match(labels, v2$AREA),]
        write.csv(h2, paste('../CSV', inst$folder, 'tables', inst$common, k, paste0('h_',rc[j],'.csv'), sep = '/'), row.names = F)
        write.csv(v2, paste('../CSV', inst$folder, 'tables', inst$common, k, paste0('v_',rc[j],'.csv'), sep = '/'), row.names = F)
      }
      max.membership <- memberships %>% dplyr::group_by(AREA) %>% dplyr::summarise(COMMSHIP = COMMSHIP[which(P == max(P))], P = max(P))
      max.membership <- max.membership[match(labels, max.membership$AREA),]
      write.csv(max.membership, paste('../CSV', inst$folder, 'tables', inst$common, k, "max.csv", sep = '/'), row.names = F)
    } 
  } else
    print("No creating tables")
}

save.tables.wsbm <- function(K, inst, labels, memberships, on=T){
  if (on){
    print("Saving WSBM table")
    print("----> Warning: Choose the right partition")
    library(magrittr)
    for (k in K){
      memberships <- memberships %>% unlist() %>% as.numeric()
      data.frame(AREA=labels, COMMSHIP=memberships, P=1, DIR="WSBM") %>% write.csv(paste('../CSV', inst$folder, 'tables', inst$common, k, paste(inst$suffix,"csv", sep = "."), sep = '/'))
    }
  
  } else{
    print("No WSBM table")
  }
}

make.hier <- function(K, net, net.cluster, regions, labels, membership, inst, on=T){
  if (on){
    print("Printing Hierarchical Edge Bundle")
    source("functions/plot_hier_edg_bundle.R")
    for (k in K){
      plot.hier.edg.bundle(k, net, net.cluster, regions, labels, membership, path=inst$plot, foldername=inst$mainfolder, subfolder="HEB", filename=paste(k, inst$suffix))   
    }
  } else{
    print("No hierarchical edge bundle")
  }
}

make.bar.block.sim <- function(R, K, net, net.cluster, inst, on=T){
  if (on){
    print("Plotting bar-block-sim-v1")
    source("functions/plot_bar_block_sim.R")
    for (r in R){
      for (k in K){
        plot.bar.block.sim(r, k, net, net.cluster, path=inst$plot, foldername=inst$mainfolder, subfolder="BlockSim", filename=sprintf("k_%i_r_%i", k, r))
      }
    }
  } else
    print("No bar-block-sim_v1")

}

make.densities <- function(net, labels, membership, on=T){
  if (on){
    print("Plotting block densities")
    source("functions/plot_densities.R")
    plot.densities(net, labels, membership)
  } else
    print("No densities")
}

make.Astats.commship <- function(net, hcluster, k, labels, memberships, inst, on=T){
  if (on){
    source("functions/plot_stats_commship.R")
    print("* Plotting different statistical plots")
    plot.area.commship(net, hcluster, k, labels, memberships, path = inst$plot, foldername = inst$mainfolder, subfolder = "stats", filename = paste("area_size","commship",k, inst$suffix, sep = "_"))
  } else
    print("No Astat commship")
}

make.Rstats.commship <- function(net, hcluster, k, labels, memberships, inst, on=T){
  if (on){
    source("functions/plot_stats_commship.R")
    print("* Creating different statistical plots")
    plot.region.commship(net, hcluster, k, labels, memberships, path = inst$plot, foldername = inst$mainfolder, subfolder = "stats", filename = paste("region_size","commship",k, inst$suffix, sep = "_"))
  } else
    print("No Rstat commship")
}

save.lincomm.matrix <- function(net, hcluster, k, inst, save=T){
  if (save){
    print("* Saving Commship matrix in WSBM folder")
    source("functions/format_lincomm.R")
    net <- format.lincomm(net, hcluster, k)
    source("functions/df_to_adj.R")
    net <- with(net, data.frame(source=source, target=target, weight=commship))
    net <- net %>% df.to.adj()
    write.csv(net, paste("../WSBM/CD/CSV/commships", inst$folder, inst$common,"4_commships_sp.csv", sep = "/"), row.names = F)
  } else
    print("Not saving link community matrix")
}

make.wdiscommship <- function(net, hcluster, k, inst, on=T){
  if (on){
    source("functions/plot_stats_commship.R")
    print("* Creating different statistical plots")
    plot.wdis.commship(net, hcluster, k, path=inst$plot, foldername=inst$mainfolder, subfolder="stats", filename="wdis_%s" %>% sprintf(k))
  } else
    print("No wdis commship")
}

main <- function(inst){
  
  library(magrittr)
  
  linkage <- 'average'
  index <- 'jaccp'
  EC <- F
  mode <- 'ALPHA'
  nlog10 <- T
  tag <- paste(inst$model, inst$distances, inst$folder, sep = "_")  ## zz_model_NONULL_tracto2016_merged   normal_tracto2016_merged
  fln.name <- 'fln'
  
  source('functions/model_name_analysis.R')
  filename <- model.name(index, linkage, mode, nlog10, EC, tag)
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
  coords <- netx$fln.coords
  regions <- netx$regions
  
  if (nlog10){
    net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
  }
  
  ### Load WSBM partition
  print("Warning: Be careful choosing the right partition")
  node.membership <- read.csv(inst$wsbm, header = F) %>% unlist() %>% as.numeric()
  source("functions/arrange_memberships_values.R")
  node.membership <- node.membership %>% arrange.memberships.values()
  
  print("** Arrange memberships with respect to a reference")
  source("functions/arrange_memberships_reference.R")
  node.membership <- arrange.memberships.reference(node.membership, inst$R)
  
  # print("** Arrange membership with some manual choices")
  # source("functions/arrange_membership_manual.R")
  # node.membership <- arrange.membership.manual(node.membership)
  
  ### Load data
  net.cluster <- readRDS(sprintf('../RDS/%s/hclust/hierarchical_clustering.rds', filename))
  hclust.features <- readRDS(sprintf('../RDS/%s/dndrgm/single_linkage_features.rds', filename))
  # ham.merde <- load.merde(hclust.features, nodes, labels, net, net.cluster, leaves, filename, save=T)
  
  ### Plots start here:
  make.parameters(inst, hclust.features, on=F)
  make.cascade(net, ham.merde, nodes, labels, regions, inst, on=F)
  make.mountain(ham.merde, labels, nodes, coords, regions, inst, on=F)
  make.circular(4, hclust.features, ham.merde, inst, regions, on=F)
  make.dend.matrix(Ns, features, ham.merde, net.cluster, nodes, inst, regions, on=F)
  make.link.comm(4, net, nodes, labels, net.cluster, leaves, regions, node.membership, inst, on=F)
  make.bars(107:2, net, ham.merde, net.cluster, labels, regions, inst, on=F)
  make.hier(4, net, net.cluster, regions, labels, node.membership, inst, on=F)
  make.bar.block.sim(6, 4, net, net.cluster, inst, on=F)
  make.densities(net, labels, node.membership, on=F)
  make.Astats.commship(net, net.cluster, 4, labels, node.membership, inst, on=F)
  make.Rstats.commship(net, net.cluster, 4, labels, node.membership, inst, on=F)
  make.wdiscommship(net, net.cluster, 4, inst, on=T)
  
  ## Table part:
  save.tables(34, net, net.cluster, labels, inst, on=F)
  save.tables.wsbm(4, inst, labels, node.membership, on=F)
  save.lincomm.matrix(net, net.cluster, 4, inst, save=F)

}

### HEAD ###

data <- "commships"
folder <- 'merged'
distances <- 'tracto2016'
model <- 'zz_model'
csv.name <- 'fln'
MAN <- "MAX"
R <- 6
labels.name <- 'imputation_labels'
wsbm.name <- paste("al_0/K",MAN, sprintf("k_%i.csv",R), sep = "/")
suffix <- paste("wsbm_al_0",MAN,R, sep = "_")


common.path <- paste(distances, model, sep = '/')
csv.path <- paste(folder, 'imputation', common.path, paste0(csv.name,'.csv'), sep = '/') # merged/imputation/tracto2016/zz_model/   91x40
labels.path <- paste(folder, 'labels', common.path, paste0(labels.name,'.csv'), sep = '/')
plot.path <- '../plots'
wsbm.path <- paste("../WSBM/CD/CSV/labels",data,folder,common.path, "NOSP", wsbm.name, sep = "/")

path.list <- list(csv=csv.path, 
                  labels=labels.path, 
                  plot=plot.path, 
                  common=common.path,
                  folder=folder,
                  model=model,
                  distances=distances,
                  wsbm=wsbm.path,
                  R=R,
                  suffix=suffix)

main(path.list)

#### Create function for densities; played a bit with save.tables.wsbm; create Rstats and Astats, making the wsbm partition unique to all functions;
#### 


### Take rhamreas between different edge communities
# rhareas_R1_B1 <- read.csv("../CSV/merged/geometry/tracto2016/zz_model/areas_R1_B1.csv", header = F) %>% unlist() %>% as.numeric()
# rhareas_R1_B2 <- read.csv("../CSV/merged/geometry/tracto2016/zz_model/areas_R1_B2.csv", header = F) %>% unlist() %>% as.numeric()
# rhareas_R2_B2 <- read.csv("../CSV/merged/geometry/tracto2016/zz_model/areas_R2_B2.csv", header = F) %>% unlist() %>% as.numeric()
# areas <- data.frame(a=c(), b=c())
# areas <- rbind(areas, data.frame(a=rhareas_R1_B1, b="R1_B1"))
# areas <- rbind(areas, data.frame(a=rhareas_R1_B2, b="R1_B2"))
# areas <- rbind(areas, data.frame(a=rhareas_R2_B2, b="R2_B2"))
# p <- ggplot2::ggplot(areas, ggplot2::aes(x=a, fill=b))+
#   ggplot2::geom_histogram(ggplot2::aes(y=..density..), bins = 30, color="gray", alpha=0.5)
# png("../plots/merged/tracto2016/zz_model/AVERAGE_full_l10/geometry/areas.png", width = 8, height = 5, units = "in", res = 200)
# print(p)
# dev.off()
