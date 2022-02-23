dis.to.2d <- function(path, dis, N){
  
  source(sprintf('%s/R/functions/adj_to_df.R', path))
  
  dis <- adj.to.df(dis)
  dis <- dis[dis$target <= N & dis$source <= N,]
  dis <- dis[dis$source != dis$target,]
  dis <- dis[dis$source < dis$target,]
  
  dis.graph <- igraph::graph_from_data_frame(dis[,1:2], directed = F)
  dis.coords <- igraph::layout_with_kk(dis.graph, weights = dis$weight)
  
  write.csv(dis.coords, sprintf('%s/CSV/dis_2d.csv', path), row.names = F)
  
}
