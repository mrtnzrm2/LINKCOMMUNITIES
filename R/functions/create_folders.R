create.folders <- function(filename){
  dir.create(sprintf('../RDS/%s/', filename), showWarnings = F)
  dir.create(sprintf('../RDS/%s/similarity/', filename), showWarnings = F)
  dir.create(sprintf('../RDS/%s/dndrgm/', filename), showWarnings = F)
  dir.create(sprintf('../RDS/%s/hclust/', filename), showWarnings = F)
}
