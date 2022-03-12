get.infra <- function(fln.labels){
  infra  <- read.csv("../CSV/Granular/infra.csv")
  source("functions/format_label_x.R")
  col.labs <- colnames(infra ) %>% format.label.x()
  col.labs <- col.labs[2:52] %>% tolower()
  source("functions/format_labels.R")
  row.labs <- infra [,1] %>% format.labels() %>% tolower()
  infra  <- infra [,2:52]
  colnames(infra ) <- col.labs
  rownames(infra ) <- row.labs
  infra  <- infra [fln.labels, fln.labels[1:50]] %>% as.matrix()
  return(infra )
}