#' FLN to 2d coordinates
#' 
#'Need fln to be ordered with N first elements to be the areas from the edge-complete part.
#'
#' @param path String: path of the main folder
#' @param fln MxN matrix: Matrix with the fln information
#' @param N Integer: length of the edge-complete graph
#'
#' @return Void
#' @export
#'
#' @examples
#' 
fln.to.2d <- function(path, fln, N){
  
  source(sprintf('%s/R/functions/adj_to_df.R', path))

  fln <- adj.to.df(fln)
  fln <- fln[fln$target <= N & fln$source <= N,]
  fln <- fln[!is.na(fln$weight),]
  
  fln.graph <- igraph::graph_from_data_frame(fln[,1:2], directed = T)
  fln.coords <- igraph::layout_with_kk(fln.graph, weights = -log10(fln$weight))
  fln.coords <- igraph::layout_with_fr(fln.graph, weights = -1/(fln$weight), coords = fln.coords)
  
  write.csv(fln.coords, sprintf('%s/CSV/fln_2d.csv', path), row.names = F)
  
}
