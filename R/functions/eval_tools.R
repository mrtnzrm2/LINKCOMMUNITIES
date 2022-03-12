norm.pred.train <- function(pred, mat, id){
  source("functions/adj_to_df.R")
  Ns <- dim(mat)
  mat <- mat %>% adj.to.df()
  nonzero <- !mat$weight == 0
  pred.mat <- mat
  pred.mat$weight[pred.mat$weight != 0] <- pred
  source("functions/df_to_adj.R")
  pred.mat <- pred.mat %>% df.to.adj()
  pred.mat[pred.mat == 0] <- NA
  pred.mat <- 10^(pred.mat-7)
  pred.mat[is.na(pred.mat)] <- 0
  pred.mat <- pred.mat/col.sum(pred.mat, Ns[2])
  pred.mat[pred.mat==0] <- NA
  pred.mat <- log10(pred.mat) + 7
  pred.mat[is.na(pred.mat)] <- 0
  pred.mat <- pred.mat %>% adj.to.df()
  pred.mat <- pred.mat$weight[nonzero]
  return(pred.mat)
}

norm.pred.test <- function(pred, mat, id){
  Ns <- dim(mat)
  pred.mat <- matrix(0, Ns[1], Ns[2])
  for (j in 1:Ns[2]){
    for (i in 1:Ns[1]){
      pred.mat[i,j] <- pred[(j-1)*Ns[1] + i]
    }
  }
  pred.mat[pred.mat == 0] <- NA
  pred.mat <- 10^(pred.mat-7)
  pred.mat[is.na(pred.mat)] <- 0
  pred.mat <- pred.mat/col.sum(pred.mat, Ns[2])
  pred.mat[pred.mat==0] <- NA
  pred.mat <- log10(pred.mat) + 7
  pred.mat[is.na(pred.mat)] <- 0
  source("functions/adj_to_df.R")
  mat <- mat %>% adj.to.df()
  nonzero <- !mat$weight == 0
  pred.mat <- pred.mat %>% adj.to.df()
  pred.mat <- pred.mat$weight[nonzero]
  return(pred.mat)
}

get.nonzero <- function(w, mat){
  source("functions/adj_to_df.R")
  mat <- mat %>% adj.to.df()
  nonzero <- !mat$weight == 0
  return(w[nonzero])
}