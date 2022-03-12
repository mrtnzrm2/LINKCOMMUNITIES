plot.performance <- function(net, aik.supra, aik.infra, aik, supra, infra, labels, path="", foldername="", subfolder="", on=T){
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  source("functions/get_tracto2016.R")
  distance <- get.tracto2016(labels)
  source("functions/adj_to_df.R")
  distance <- distance %>% adj.to.df()
  distance <- distance[distance$target <= 50,]
  
  source("functions/df_to_adj.R")
  net <- net %>% df.to.adj()
  net <- net %>% adj.to.df()
  net <- net[net$target <= 50,]
  
  zeros <- net$weight == 0
  aik.supra <- aik.supra[aik.supra$target <= 50,]
  aik.infra <- aik.infra[aik.infra$target <= 50,]
  aik <- aik[aik$target <= 50,]
  
  supra <- supra$weight[!zeros]
  infra <- infra$weight[!zeros]
  distance <- distance$weight[!zeros]
  w <- net$weight[!zeros]
  aik <- aik$weight[!zeros]
  
  w <- w[-1290]
  supra <- supra[-1290]
  infra <- infra[-1290]
  distance <- distance[-1290]
  aik <- aik[-1290]
  
  data <- data.frame(supra=supra-mean(supra), infra=infra-mean(infra), distance=distance-mean(distance), w = w, aik=aik-mean(aik))
  model <- lm(w ~., data = data)
  
  p <- performance::check_model(model)
  
  png(sprintf("%s/%s/%s/check_model.png", path, foldername, subfolder), width = 10, height = 10, units = 'in', res = 200)
  print(p)
  dev.off()

}