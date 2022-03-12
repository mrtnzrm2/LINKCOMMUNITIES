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

red.laboratory <- function(A, memberships, inst, on=T){
  ### RED laboratory
  if (on){
    print("Red laboratory working")
    red <- A
    red$source.commship <- memberships[red$source]
    red$target.commship <- memberships[red$target]
    src <- 1:6
    tgt <- 1:6
    com <- 1
    red$flag <- ifelse(red$source.commship %in% src & red$target.commship %in% tgt & red$weight %in% com, 1, 0)
    red <- with(red, data.frame(source=source, target=target, weight=flag))
    red <- df.to.adj(red)
    p <- plot.matrix(red)
    
    png(paste(inst$plot, inst$mainfolder, "RnB", paste0("red_src_",paste0(src, collapse = ""),"_tgt_", paste0(tgt, collapse = ""), "_lcom_", paste0(com, collapse = ""), ".png"), sep = "/"), width = 6, height = 5, units = "in", res = 200)
    print(p)
    dev.off()
    
    write.csv(red, paste("../CSV", inst$folder, "RnB", inst$common, paste0("red_src_", paste0(src, collapse = ""), "_tgt_", paste0(tgt, collapse = ""), "_lcom_", paste0(com, collapse = ""),".csv"), sep = "/"), row.names = F)
  } else
    print("Red laboratory day-off!")
  
}

blue.laboratory <- function(A, memberships, inst, on=T){
  ### RED laboratory
  if (on){
    print("Blue laboratory working")
    red <- A
    red$source.commship <- memberships[red$source]
    red$target.commship <- memberships[red$target]
    
    src <- 1:6
    con.src <- -1
    tgt <- 1:6
    con.tgt <- -1
    com <- 1
    flag.1 <- ifelse(red$source.commship %in% src & red$target.commship %in% tgt & red$weight %in% com, 1, 0)
    flag.1 <- flag.1 * ifelse(red$source.commship %in% con.src & red$target.commship %in% con.tgt & red$weight %in% com, 0, 1)
    
    
    # src <- 1:6
    # con.src <- -1
    # tgt <- 1:6
    # con.tgt <- -1
    # com <- 2
    # flag.2 <- ifelse(red$source.commship %in% src & red$target.commship %in% tgt & red$weight %in% com, 1, 0)
    # flag.2 <- flag.2 * ifelse(red$source.commship %in% con.src & red$target.commship %in% con.tgt & red$weight %in% com, 0, 1)
    
    red$flag <- flag.1 #+ flag.2
    red <- with(red, data.frame(source=source, target=target, weight=flag))
    red <- df.to.adj(red)
    p <- plot.matrix(red)
    
    suffix <- paste0(com, collapse = "")
    if (!-1 %in% con.src)
      suffix <- paste(suffix, con.src %>% paste0(collapse = ""), sep = "0")
    if (!-1 %in% con.tgt)
      suffix <- paste(suffix, con.tgt %>% paste0(collapse = ""), sep = "0")
    
    suffix <- 'special'
    
    png(paste(inst$plot, inst$mainfolder, "RnB", paste0("blue_src_", paste0(src, collapse = ""),"_tgt_", paste0(tgt, collapse = ""), "_lcom_", suffix, ".png"), sep = "/"), width = 6, height = 5, units = "in", res = 200)
    print(p)
    dev.off()
    
    write.csv(red, paste("../CSV", inst$folder, "RnB", inst$common, paste0("blue_src_", paste0(src, collapse = ""), "_tgt_", paste0(tgt, collapse = ""), "_lcom_", suffix, ".csv"), sep = "/"), row.names = F)
  } else
    print("Blue laboratory day-off!")
  
}

arrange.membership.manual <- function(id){
  # 3 -> 1 1 -> 3 5 -> 6 6 -> 5
  for (i in 1:length(id)){
    if (id[i] == 1)
      id[i] <- 3
    else if (id[i] == 3)
      id[i] <- 1
    else if (id[i] == 5)
      id[i] <- 6
    else if (id[i] == 6)
      id[i] <- 5
  }
  
  return(id)
}

main <- function(inst, k){
  library(magrittr)
  source("functions/adj_to_df.R")
  source("functions/df_to_adj.R")
  source("functions/plot_matrix.R")
  
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
  print(inst$mainfolder)
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
  
  if (nlog10)
    net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
  
  ### Load hclust from selected dataset
  net.cluster <- readRDS(sprintf('../RDS/%s/hclust/hierarchical_clustering.rds', filename))
  
  ### In A, each edge has its link commmunity membership
  A <- commship.n.filter(net, net.cluster, k)
  A <- with(A, data.frame(source=source, target=target, weight=commship))
  A <- df.to.adj(A)
  
  print("Warning: Be careful choosing the right partition")
  node.membership <- read.csv("../WSBM/CD/CSV/labels/merged/tracto2016/zz_model/4_r_6_3.csv", header = F) %>% as.matrix()
  node.membership <- arrange.membership.manual(node.membership)
  ### Reorder using node memberships
  A <- A[order(node.membership), order(node.membership)]
  node.membership <- node.membership[order(node.membership)]
  A <- A  %>% adj.to.df()
  A <- A[A$weight > 0,]
  
  ### Go to lab!!
  red.laboratory(A, node.membership, inst, on = T)
  blue.laboratory(A, node.membership, inst, on = F)
  
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
                  folder=folder,
                  model=model,
                  distances=distances)

main(path.list, 4)