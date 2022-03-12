compute.aki.aik <- function(net, nodes, mode, nt, inst, suffix="", save=T){
  if (save){
    ### Compute AKI and AIK
    source("functions/aki.R")
    source("functions/aik.R")
    source("functions/jaccard_p_fast.R")
    
    ns <- 107
    
    net <- net[net$target <= nt,]
    
    aki.matrix <- aki(net, mode=mode)
    aik.matrix <- aik(net, mode=mode)
    
    AIK <- matrix(0, nrow = ns, ncol = ns)
    AKI <- matrix(0, nrow = nt, ncol = nt)
    
    for (i in 1:ns){
      for (j in 1:ns){
        if (i < j){
          AIK[i,j] <- jaccard.p.fast(aik.matrix[i,], aik.matrix[j,])

        }
      }
    }
    
    for (i in 1:nt){
      for (j in 1:nt){
        if (i < j){
          AKI[i,j] <- jaccard.p.fast(aki.matrix[i,], aki.matrix[j,])
        }
      }
    }
    
    AIK <- AIK + t(AIK)
    AKI <- AKI + t(AKI)
    AKI[diag(nt)==1] <- -1
    AIK[diag(ns)==1] <- -1
    
    write.csv(AKI, paste("../CSV", inst$folder, "similarity", inst$common, "aki_compact%s.csv" %>% sprintf(suffix), sep = "/"), row.names = F)
    write.csv(AIK, paste("../CSV", inst$folder, "similarity", inst$common, "aik_compact%s.csv" %>% sprintf(suffix), sep = "/"), row.names = F)
  } else{
    ### Load AKI and AIK
    AKI <- read.csv(paste("../CSV", inst$folder, "similarity", inst$common, "aki_compact%s.csv" %>% sprintf(suffix), sep = "/")) %>% as.matrix()
    AIK <- read.csv(paste("../CSV", inst$folder, "similarity", inst$common, "aik_compact%s.csv" %>% sprintf(suffix), sep = "/")) %>% as.matrix()
  }
  A <- list(AKI=AKI, AIK=AIK)
  return(A)
}