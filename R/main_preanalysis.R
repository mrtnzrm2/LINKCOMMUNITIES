make.aks.matrix <- function(AKS, bourder, colors, labels, regions, inst, preffix="", suffix="", on=T){
  if (on){
    print("** Plotting Aks matrix")
    source("functions/plot_aks_matrix.R")
    plot.aks.matrix(AKS, bourder, colors, labels, regions, path=inst$plot, foldername=inst$mainfolder, subfolder="Preprocessing", preffix=preffix, suffix=suffix)
  } else
    print("** No Aks matrix")
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

make.fln.histogram <- function(net, inst, on=T){
  if (on){
    print("Plotting fln histogram")
    source("functions/plot_fln_histogram.R")
    plot.fln.histogram(net, path=inst$plot, foldername=inst$mainfolder, subfolder="Preprocessing")
  } else
    print("No fln histogram")
}

make.performance <- function(df, labels, inst, on=T){
  if (on){
    print("** Plotting performance")
    source("functions/plot_performance.R")
    plot.performance(df$net, df$aik.supra, df$aik.infra, df$aik, df$supra, df$infra, labels, path=inst$plot, foldername=inst$mainfolder, subfolder="Preprocessing")
  } else
    print("**No performance")
}

col.sum <- function(A, N){
  x <- c()
  for ( i in 1:N){
    x <- c(x, sum(A[,i], na.rm = T))
  }
  return(x)
}

get.datasets <-  function(net, labels, nt, inst, filename=""){
  source("functions/adj_to_df.R")
  source("functions/df_to_adj.R")
  # Ensemble data ----
  net <- net[net$target <= nt,] %>%
    df.to.adj() %>% 
    adj.to.df() %>% 
    dplyr::mutate(weight=ifelse(weight > 0, weight, NA))
  source("functions/get_supra.R")
  supra <- get.supra(labels)
  supra[supra == 0] <- NA
  # supra <- supra/col.sum(supra, nt)
  supra <- adj.to.df(supra) %>% 
    dplyr::mutate(weight=log10(weight))
  source("functions/get_infra.R")
  infra <- get.infra(labels) 
  infra[infra==0] <- NA
  # infra <- infra/col.sum(infra, nt)
  infra <- adj.to.df(infra) %>% 
    dplyr::mutate(weight=log10(weight))
  source("functions/get_tracto2016.R")
  distance <- get.tracto2016(labels) 
  distance <- distance[,1:nt] %>%
    adj.to.df()

  data <- dplyr::tibble(
    w=net %>%
      dplyr::pull(weight),
    supra=supra %>%
      dplyr::pull(weight),
    infra=infra %>%
      dplyr::pull(weight),
    distance=distance %>%
      dplyr::pull(weight)
    )
  
  data <- do.call(data.frame, data %>% lapply(function(x) replace(x, x == 0, NA)))
  
  print(skimr::skim(data))
  
  p <- data %>% dplyr::select_if(is.numeric) %>%
    GGally::ggpairs()
  
  png("%s/%s/PreAnalysis/describe_%s.png" %>% sprintf(inst$plot, inst$folder, filename), width = 8, height = 8, res = 200, units = "in")
  print(p)
  dev.off()
  
  return(list(net=net, supra=supra, infra=infra, distance=distance))
}

get.mean.similarity <- function(df, ntrn, ntgt, nsrc=107){
  
  net <- df$net
  AIK <- df$aik
  AKI <- df$aki
  AIK[diag(nsrc) == 1] <- NA
  AKI[diag(ntrn) == 1] <- NA
  similarity <- matrix(0, nrow = nsrc, ncol = ntgt)
  
  for (i in 1:ntrn){
    for (j in 1:ntrn){
      if (i != j){
        jacp.aik <- AIK[i,j]
        jacp.aki <- AKI[i,j]
        similarity[i,j] <- mean(c(jacp.aki, jacp.aik))
      }
    }
  }
  
  for (i in (ntrn+1):nsrc){
    for (j in 1:ntrn){
      Ni <- which(net[i,] > 0)
      Nj <- which(net[j,] > 0)
      NiNj <- Ni[Ni %in% Nj]
      jacp.aik <- AIK[i,j]
      jacp.aki <- mean(AKI[j, NiNj], na.rm = T)
      similarity[i,j] <- mean(c(jacp.aki, jacp.aik))
    }
  }
  return(similarity)
}

transform.datasets <- function(dataset, nodes, nt, inst, filename=""){
  # Define datasets ----
  net <- dataset$net  %>% df.to.adj()
  net[is.na(net)] <- 0
  supra <- dataset$supra  %>% df.to.adj()
  supra[is.na(supra)] <- 0
  infra <- dataset$infra  %>% df.to.adj()
  infra[is.na(infra)] <- 0
  distance <- dataset$distance
  # Compute Aks ----
  source("functions/compute_aik.R")
  source("functions/compute_aki.R")
  net.aik <- compute.aik(net %>% adj.to.df(), nodes)
  net.aki <- compute.aki(net %>% adj.to.df(), nt)
  supra.aik <- compute.aik(supra %>% adj.to.df(), nodes)
  supra.aki <- compute.aki(supra %>% adj.to.df(), nt)
  infra.aik <- compute.aik(infra %>% adj.to.df(), nodes)
  infra.aki <- compute.aki(infra %>% adj.to.df(), nt)
  # Compute similarity matrices ----
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  net.sim <- get.mean.similarity(
    list(
      net=net,
      aik=net.aik, 
      aki=net.aki
      ), nt, nt) %>% adj.to.df()
  supra.sim <- get.mean.similarity(
    list(
      net=supra, 
      aik=supra.aik,
      aki=supra.aki
      ), nt, nt) %>% adj.to.df()
  infra.sim <- get.mean.similarity(
    list(
      net=infra,
      aik=infra.aik,
      aki=infra.aki
      ), nt, nt) %>% adj.to.df()
  
  data <- dplyr::tibble(
    w=net %>% adj.to.df() %>% dplyr::pull(weight),
    supra=supra.sim %>% dplyr::pull(weight),
    infra=infra.sim %>% dplyr::pull(weight),
    sim=net.sim %>% dplyr::pull(weight),
    distance=distance %>% dplyr::pull(weight)
  )
  
  data <- do.call(data.frame, data %>% lapply(function(x) replace(x, x == 0, NA)))
  
  print(skimr::skim(data))
  
  p <- data %>% dplyr::select_if(is.numeric) %>%
    GGally::ggpairs()
  
  png("%s/%s/PreAnalysis/describe_%s.png" %>% sprintf(inst$plot, inst$folder, filename), width = 8, height = 8, res = 200, units = "in")
  print(p)
  dev.off()
}

main <- function(inst){
  library(magrittr)
  source("functions/adj_to_df.R")
  
  mode <- 'ALPHA'
  nlog10 <- T
  fln.name <- 'fln'
  nt <- 50
  nodes <- 107
  suffix <- "_EC"
  
  inst$mainfolder <- paste(inst$folder, "PreAnalysis", sep = "/")
  dir.create(sprintf('%s/%s', inst$plot, inst$mainfolder), showWarnings = F)
  
  # Load data ----
  source('functions/load_net.R')
  netx <- load.net(inst)
  net <- netx$net
  leaves <- netx$leaves
  labels <- netx$labels
  source("functions/format_labels.R")
  labels <- labels %>% format.labels()
  
  if (nlog10){
    net$weight[net$weight != 0] <- log10(net$weight[net$weight != 0]) + 7
  }
  
  source("functions/get_regions_colors.R")
  cr <- get.regions.colors(labels, nt)
  
  # Analyzing data ----
  datasets <- get.datasets(net, labels, nt, inst, filename="default_data")
  transform.datasets(datasets, nodes, nt, inst, filename = "transform_data")
  # Plotting ----

}


### HEAD ###
folder <- 'merged'
distances <- 'original'
model <- 'normal'
csv.name <- 'fln'
labels.name <- 'rlabels'  # imputation_labels  rlabels

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