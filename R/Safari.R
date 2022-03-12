commship.n.filter <- function(net, hclust, k){
  source('functions/assign_commship.R')
  net <- assign.commship(net, hclust, k, 0)
  source('functions/linkcomm_parameters.R')
  est.para <- linkcomm.parameters(net)
  net$commship[which(net$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net$commship[which(net$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  net <- net[net$commship > 0,]
  return(net)
}

make.link.comm <- function(Ns, net, nodes, labels, cluster, leaves, regions, inst, on=T){
  if (on){
    print("Plotting link communities matrices")
    source('functions/plot_lincomm_matrix.R')
    for (k.2 in Ns){
      plot.lincomm.matrix(net[net$source <= nodes & net$target <= nodes,], cluster,
                          leaves, nodes, k.2, 0,
                          labels, regions, plt=T, subfolder='Matrices',
                          path=inst$plot,
                          filename = sprintf('%s', k.2),
                          foldername = inst$mainfolder)
    } 
  } else
    print("No link communities matrix")
}

load.merde <- function(process, nodes, labels, net, net.cluster, leaves, filename, save=T){
  if (save){
    print("Loading the merde")
    dir.create(sprintf('../RDS/%s/merde/', filename), showWarnings = F)
    source('functions/la_arbre_a_merde.R')
    ham.merde <- la.arbre.a.merde(process, nodes, labels[1:nodes], net, net.cluster, leaves)
    # saveRDS(ham.merde,  paste('../RDS', filename, 'merde', 'node_agglomeration.rds', sep = '/'))
  } else{
    print("Not loading the merde")
    ham.merde <- readRDS(paste('../RDS', filename, 'merde', 'node_agglomeration.rds', sep = '/'))
  }
  return(ham.merde)
}

library(magrittr)



# nmi.vector <- rep(0, 14)
# 
# for (i in 2:15){
#   wsbm.zz_model <- read.csv("../WSBM/CD/CSV/labels/commships/merged/tracto2016/zz_model/al_0/K/MAX/k_%i.csv" %>% sprintf(i)) %>% unlist() %>% as.numeric()
#   wsbm.normal <- read.csv("../WSBM/CD/CSV/labels/commships/merged/tracto2016/normal/al_0/K/MAX/k_%i.csv" %>% sprintf(i)) %>% unlist() %>% as.numeric()
#   nmi.vector[i-1] <- aricode::NMI(wsbm.zz_model, wsbm.normal, variant = "joint")
# }
# 
# print(nmi.vector)


# path <- '/Users/jmarti53/Documents/ND/LINKCOMMUNITIES/MAC/220126'
# path.save <- paste(path, 'plots', sep = '/')
# foldername <- ''
# 
# source('functions/df_to_adj.R')
# source("functions/adj_to_df.R")
# source('functions/load_net.R')
# source("functions/plot_matrix.R")
# source('functions/jaccard_p_fast.R')
# 
# library(magrittr)
# 
# folder <- 'merged'
# distances <- 'tracto2016'
# model <- 'zz_model'
# csv.name <- 'fln_4_r_6_3'  # fln_4_r_6_3  fln
# labels.name <- 'imputation_labels'
# 
# common.path <- paste(distances, model, sep = '/')
# csv.path <- paste(folder, 'imputation', common.path, paste0(csv.name,'.csv'), sep = '/') # merged/imputation/tracto2016/zz_model/   91x40
# labels.path <- paste(folder, 'labels', common.path, paste0(labels.name,'.csv'), sep = '/')
# plot.path <- '../plots'
# 
# inst <- list(csv=csv.path, 
#             labels=labels.path, 
#             plot=plot.path, 
#             common=common.path,
#             folder=folder,
#             distances=distances,
#             model=model)
# 
# linkage <- 'average'
# index <- 'jaccp'
# EC <- F
# mode <- 'ALPHA'
# nlog10 <- T
# tag <- paste(inst$model, inst$distances, inst$folder, sep = "_")  ## zz_model_NONULL_tracto2016_merged   normal_tracto2016_merged
# 
# source('functions/model_name_analysis.R')
# filename <- model.name(index, linkage, mode, nlog10, EC, tag)
# print(filename)
# source('functions/sformat.R')
# inst$mainfolder <- paste(inst$folder, inst$common, paste(toupper(linkage), 'full','l10', sep = "_"), sep = "/")
# dir.create(sprintf('%s/%s', inst$plot, inst$mainfolder), showWarnings = F)
# 
# #### Load network
# source('functions/load_net.R')
# netx <- load.net(inst)
# net <- netx$net
# nodes <- netx$nodes
# leaves <- netx$leaves
# labels <- netx$labels
# coords <- netx$fln.coords
# regions <- netx$regions
# 
# if (nlog10){
#   net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
# }
# 
# ### Load data
# net.sim <- read.csv(paste('../CSV', inst$folder,'similarity', inst$common,'sim_full_l10_4_r_6_3.csv', sep = "/"), header = F) %>% as.matrix() %>% t()  #   sim_full_l10_4_r_6_3.csv
# net.sim[net.sim == -1] <- 0
# net.sim <- net.sim + t(net.sim)
# net.sim[net.sim == 0] <- NA
# 
# net.dis <- -log10(net.sim)
# net.dis[is.na(net.dis)] <- max(net.dis, na.rm = T) + 1
# 
# ### Load hclust from selected dataset
# net.cluster <- hclust(as.dist(net.dis), method = linkage)
# 
# ### In A, each edge has its link commmunity membership
# nat <- commship.n.filter(net, net.cluster, 4)
# nat <- with(nat, data.frame(source=source, target=target, weight=commship))
# nat$id <- 1:nrow(nat)
# list.1 <- c(1,1,1,1,2,2,2,3,3,4)
# list.2 <- c(1,2,3,4,2,3,4,3,4,4)
# for (i in 1:length(list.1)){
#   id.1 <- nat$id[nat$weight==list.1[i]]
#   id.2 <- nat$id[nat$weight==list.2[i]]
#   m <- mean(net.dis[id.1, id.2])
#   print(m)
# }
# 
# make.link.comm(4, net, nodes, labels, net.cluster, leaves, regions, inst, on=T)
