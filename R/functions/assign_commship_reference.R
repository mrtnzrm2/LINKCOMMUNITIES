assign.commship.reference <- function(net, net.cluster, k) {
  ref.hclust <- readRDS("../RDS/jaccp_average_ALPHA_normal_tracto2016_merged_l10/hclust/hierarchical_clustering.rds")
  ref.commship <- cutree(ref.hclust, k = k)
  org.commship <- cutree(net.cluster, k = k)
  source("functions/apply_munkres.R")
  net$commship <- apply.munkres(org.commship, ref.commship, k)
  return(net)
}