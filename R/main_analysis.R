main <- function(inst){
  
  library(magrittr)

  linkage <- 'average'
  similarity.index <- 'jaccp'
  self.loop <- T
  self.loop.type <- 'ALPHA'
  simSave <- F
  estSave <- T
  EC <- F
  nlog10 <- T
  tag <- 'zz_model_NONULL_tracto2016_merged' # zz_model_NONULL_tracto2016_merged  91x40  normal_tracto2016_merged
  fln.name <- 'fln'

  source('functions/model_name_analysis.R')
  filename <- model.name(similarity.index, linkage, self.loop.type, nlog10, EC, tag)
  source('functions/create_folders.R')
  create.folders(filename)

  ### Load network
  source('functions/load_net.R')
  fln <- load.net(inst)
  net <- fln$net
  nodes <- fln$nodes
  leaves <- fln$leaves
  
  if (nlog10){
    net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
  }
  

  source('functions/sformat.R')
  if (simSave){
    source('functions/similarity_competitive.R')
    net.sim <- similarity.competitive(net, leaves, self.loop = T, mode=self.loop.type)
    saveRDS(net.sim, file.path('..', 'RDS', '{name}/similarity/link_similarity.rds' %>% sformat(list(name=filename))) )
  } else{
    # net.sim <- readRDS(file.path('..', 'RDS', '{name}/similarity/link_similarity.rds' %>% sformat(list(name=filename))) )
    net.sim <- read.csv(paste('../CSV', inst$folder,'similarity', inst$common,'sim_full_l10.csv', sep = "/"), header = F) %>% as.matrix() %>% t()
    net.sim[net.sim == -1] <- NA
  }

  net.dis <- -log10(net.sim)
  net.dis[is.na(net.dis)] <- max(net.dis, na.rm = T) + 1

  net.cluster <- hclust(as.dist(net.dis), method = linkage)
  saveRDS(net.cluster, file.path('..', 'RDS', '{name}/hclust/hierarchical_clustering.rds' %>% sformat(list(name=filename))) )

  if(estSave){
    source('functions/process_hclust.R')
    hclust.features <- process.hclust(net, net.cluster, nodes)
    saveRDS(hclust.features, file.path('..', 'RDS', '{name}/dndrgm/single_linkage_features.rds' %>% sformat(list(name=filename))) )
  } else{
    hclust.features <- readRDS(file.path('..', 'RDS', '{name}/dndrgm/single_linkage_features.rds' %>% sformat(list(name=filename))) )
  }

  print(hclust.features)
  
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

# ### FOR IMPUTATION
# get.edc <- net$source <= nodes & net$target <= nodes
# leaves <- net[get.edc, ] %>% nrow()
# net$id[get.edc] <- 1:leaves
# net$id[!get.edc] <- (leaves+1):nrow(net)
# net <- net[order(net$id),]
# ###
# 
# ### Truncate network to get the edge-complete subgraph
# if (grepl( 'EC', filename,  fixed = TRUE)){
#   net <- net[which(net$source <= max(net$target)),]
#   leaves <- nrow(net)
#   net$id  <- 1:leaves
# }
