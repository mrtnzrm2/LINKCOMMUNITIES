arrange.memberships.reference <- function(memberships, k){
  ref.membership <- read.csv("../WSBM/CD/CSV/labels/commships/merged/tracto2016/zz_model/NOSP/al_0/K/MAX/k_%i.csv" %>% sprintf(k)) %>% unlist() %>% as.numeric()
  source("functions/apply_munkres.R")
  # print(memberships)
  # print(ref.membership)
  memberships <- apply.munkres(memberships, ref.membership, k)
  # print(memberships)
  return(memberships)
}