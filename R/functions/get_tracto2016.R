get.tracto2016 <- function(fln.labels){
  tracto2016 <- read.csv("../CSV/Distances/dis_tra_11_mn_ld_2022labels.csv") %>% as.matrix() 
  tracto.labels <- tracto2016[,1] %>% tolower()
  tracto2016 <- tracto2016[,2:108] %>% as.numeric() %>% pracma::Reshape(107,107)
  tracto2016 <- tracto2016 %>% as.data.frame()
  colnames(tracto2016) <- tracto.labels
  rownames(tracto2016) <- tracto.labels
  map3d <- read.csv("../CSV/Distances/DistanceMatrix_Map3D_10feb2022_107x107.csv") %>% as.matrix()
  map3d.labels <- map3d[,1] %>% tolower()
  map3d <- map3d[,2:108] %>% as.numeric() %>% pracma::Reshape(107,107)
  map3d <- map3d %>% as.data.frame()
  colnames(map3d) <- map3d.labels
  rownames(map3d) <- map3d.labels
  
  map3d <- map3d[tracto.labels, tracto.labels] %>% as.matrix()
  
  source("functions/adj_to_df.R")
  tracto2016 <- tracto2016 %>% as.matrix() %>% adj.to.df()
  map3d <- map3d %>% adj.to.df()
  
  lm.df <- data.frame(x=map3d$weight[!is.na(tracto2016$weight)], y=tracto2016$weight[!is.na(tracto2016$weight)])
  p.lm <- lm(data = lm.df, y ~ x)
  
  source("functions/df_to_adj.R")
  tracto2016 <- tracto2016 %>% df.to.adj()
  tracto2016[diag(107) == 1] <- 0
  map3d <- map3d %>% df.to.adj()
  
  tracto.idx <- which(is.na(tracto2016), arr.ind = T)
  pred.tracto <- map3d[tracto.idx]
  
  tracto2016[tracto.idx] <- predict(p.lm, data.frame(x=pred.tracto))
  tracto2016[diag(107) == 1] <- NA
  
  tracto2016 <- tracto2016 %>% as.data.frame()
  source("functions/format_labels.R")
  tracto.labels <- format.labels(tracto.labels)
  colnames(tracto2016) <- tracto.labels
  rownames(tracto2016) <- tracto.labels
  tracto2016 <- tracto2016[fln.labels, fln.labels] %>% as.matrix()
  
  return(tracto2016)
}