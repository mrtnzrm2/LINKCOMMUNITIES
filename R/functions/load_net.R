load.net <- function(inst){
  
  library(magrittr)
  source('functions/sformat.R')
  
  net.fln <- read.csv(file.path('..', 'CSV', inst$csv))
  fln.coords<- read.csv( file.path('..', 'CSV', 'fln_2d.csv')) %>% as.matrix()
  dis.coords<- read.csv(file.path('..', 'CSV', 'dis_2d.csv')) %>% as.matrix()
  
  net.fln <- as.matrix(net.fln)
  net.fln[net.fln == 0] <- NA

  labels <- read.csv(file.path('..', 'CSV', inst$labels)) %>% as.matrix() %>% unname() %>% tolower() # 'merged', 'labels',  'imputation_zz_model_NONULL_tracto2016_clabel.csv'
  # labels <-  read.csv(sprintf('%s/CSV/CountMatrix_Summed_Main43_Periphery8_220126_arr_rlabels.csv', path)) %>% as.matrix() %>% unname() %>% tolower()
  
  leaves <- length(net.fln[!is.na(net.fln)])
  net <- matrix(0, nrow = leaves, ncol = 3)
  cnt <- 1
  ns <- nrow(net.fln)
  nt <- ncol(net.fln)
  for (i in 1:ns){
    for (j in 1:nt){
      if (!is.na(net.fln[i,j])){
        net[cnt,] <- c(i, j, net.fln[i,j])
        cnt <- cnt+ 1
      }
    }
  }

  net <- as.data.frame(net)
  colnames(net) <- c('source', 'target', 'weight')
  net$id <- 1:leaves
  nodes <- max(net$target)
  net$fln <- net$weight
  net$w <- log10(net$weight) + 7
  
  dorsal <- c( "7A", "STPc", "STPi", "STPr", "FST", "MST", "MT", "LIP" , "PIP" , "MIP" , "VIP" , "PGa", "V3A" , "DP" , "V6" , "V6A" , "V4pcLF" , "V4fpLF" ) %>% tolower()
  ventral <- c( "TH/TF", "TEav", "TEad", "TEa/mp", "TEa/ma", "TEpv", "TEpd", "TEOm", "TEO", "PERI" , "V4c" , "V4pcUF" , "V4fpUF"  ) %>% tolower()
  
  sudo.dorsal <- c('V1pcLF', 'V1fpLF') %>% tolower()
  sudo.ventral <- c('V1fpUF', 'V2fpUF', 'V2pcUF', 'V1c') %>% tolower()
  
  regions <- data.frame(AREA = labels, REGION=NA)
  regions$REGION[ labels %in% dorsal] <- 'dorsal'
  regions$REGION[ labels %in% ventral] <- 'ventral'
  regions$REGION[ labels %in% sudo.dorsal] <- 'sudo.dorsal'
  regions$REGION[ labels %in% sudo.ventral] <- 'sudo.ventral'
  regions$REGION[is.na(regions$REGION)] <- 'other'
  regions$COLOR <- NA
  regions$COLOR[regions$REGION == 'dorsal'] <- rgb(1,0.2,0,0.9)
  regions$COLOR[regions$REGION == 'ventral'] <- rgb(25,25,112, 230, maxColorValue = 255)
  regions$COLOR[regions$REGION == 'sudo.dorsal'] <- rgb(233,150,122, 230, maxColorValue = 255)
  regions$COLOR[regions$REGION == 'sudo.ventral'] <- rgb(100,149,237, 230, maxColorValue = 255)
  regions$COLOR[regions$REGION == 'other'] <- rgb(60,179,113,230, maxColorValue = 255)

  x <- list()

  x$net <- net
  x$nodes <- nodes
  # x$leaves <- length(net$source[which(net$source <= nodes)])
  x$leaves <- nrow(net)
  x$labels <- labels
  x$fln.coords <- fln.coords
  x$dis.coords <- dis.coords
  x$regions <- regions

  return(x)
}