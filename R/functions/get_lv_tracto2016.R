get.lv.tracto2016 <- function(fln.labels){
  tracto2016 <- read.csv("../CSV/Distances/LV_zz_model_tracto2016_merged.csv") %>% as.matrix() 
  tracto.labels <- read.csv("../CSV/merged/labels/original/normal/rlabels.csv") %>% as.matrix() %>% as.character() %>% tolower()
  tracto2016 <- tracto2016%>% as.numeric() %>% pracma::Reshape(107,107)
  tracto2016 <- tracto2016 %>% as.data.frame()
  source("functions/format_labels.R")
  tracto.labels <- format.labels(tracto.labels)
  colnames(tracto2016) <- tracto.labels
  rownames(tracto2016) <- tracto.labels
  tracto2016 <- tracto2016[fln.labels, fln.labels] %>% as.matrix()
  
  return(tracto2016)
}