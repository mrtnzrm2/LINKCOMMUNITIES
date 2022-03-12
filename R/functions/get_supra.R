get.supra <- function(fln.labels){
  
  supra <- read.csv("../CSV/Granular/supra.csv")
  source("functions/format_label_x.R")
  col.labs <- colnames(supra) %>% format.label.x()
  col.labs <- col.labs[2:52] %>% tolower()
  source("functions/format_labels.R")
  row.labs <- supra[,1] %>% format.labels() %>% tolower()
  supra <- supra[,2:52]
  colnames(supra) <- col.labs
  rownames(supra) <- row.labs
  supra <- supra[fln.labels, fln.labels[1:50]] %>% as.matrix()
  return(supra)
}