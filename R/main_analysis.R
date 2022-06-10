main <- function(inst) {
  library(magrittr)
  linkage <- "average"
  index <- "jaccp"
  mode <- "ALPHA"
  simSave <- F
  estSave <- T
  EC <- F
  nlog10 <- T
  fln.name <- "fln"
  tag <- paste(inst$model, inst$distances, inst$folder, sep = "_")
  source("functions/model_name_analysis.R")
  filename <- model.name(index, linkage, mode, nlog10, EC, tag)
  source("functions/create_folders.R")
  create.folders(filename)
  # Load network ----
  source("functions/load_net.R")
  fln <- load.net(inst)
  net <- fln$net
  nodes <- fln$nodes
  leaves <- fln$leaves
  ## Transform to -log10 ----
  if (nlog10)
    net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
  source("functions/sformat.R")
  # Compute similarity edge matrix ----
  if (simSave) {
    print("Calculating sim matrix using R")
    Rcpp::sourceCpp("../cpp/similarity_cooperative.cpp")
    net_sim <- similarity_cooperative(
      file.path("..", "CSV", inst$csv),
      mode,
      nodes,
      leaves
    ) %>%
     unlist() %>%
     pracma::Reshape(leaves, leaves)
    net_sim[net_sim == -1] <- NA
    print("* Finished")
    print("* Save sim matrix")
    saveRDS(
      net_sim,
      file.path(
        "..",
        "RDS",
        "{name}/similarity/link_similarity.rds" %>%
          sformat(list(name = filename))
      )
    )
    print("* Finished")
  } else{
    print("** Warning: Be careful. Load sim matrix")
    net_sim <- readRDS(
      file.path(
        "..",
        "RDS",
        "{name}/similarity/link_similarity.rds" %>%
          sformat(list(name = filename))
      )
    )
    net_sim[net_sim == -1] <- NA
    print("* Finished")
  }
  # Similarity matrix to distance matrix ----
  print("* Formatting sim matrix")
  net_dis <- -log10(net_sim)
  net_dis[is.na(net_dis)] <- max(net_dis, na.rm = T) + 1
  print("* Finished")
  print("* Getting hclust from sim matrix")
  net_cluster <- hclust(as.dist(net_dis), method = linkage)
  print("* Finished")
  print("* Save hclust")
  saveRDS(
    net_cluster,
    file.path(
      "..",
      "RDS",
      "{name}/hclust/hierarchical_clustering.rds" %>%
        sformat(list(name = filename))
    )
  )
  print("* Finished")
  source("functions/assign_commship_reference.R")
  net <- assign.commship.reference(net, net_cluster, 4)
  # Compute dendrogram heuristics ----
  print("* Analysing the hclust")
  if (estSave) {
    source("functions/fx.R")
    hclust_features <- dplyr::tibble()
    Rcpp::sourceCpp("../cpp/process_hclust.cpp")
    for (nss in c(2, 3, 4, 20)) {
      print(
        "- nss: %i" %>%
          sprintf(nss)
      )
      for (th in pracma::linspace(0.1, 0.3, 10)) {
        print(
        "--- th: %.2f" %>%
          sprintf(th)
      )
        hf <- process_hclust_fast(
          leaves, fx(net_dis),
          net$source, net$target,
          nodes, nss,
          th
        ) %>%
          t() %>%
          dplyr::as_tibble()
        colnames(hf) <- c(
          "K", "Dc",
          "PSS", "height"
        )
        hf <- hf %>%
          dplyr::mutate(
            nss = nss,
            th = th
          )
        hclust_features <- hclust_features %>%
          dplyr::bind_rows(hf)
      }
    }
    print("* Finished")
    print("* Save hclust analysis")
    saveRDS(
      hclust_features,
      file.path(
        "..",
        "RDS",
        "{name}/dndrgm/single_linkage_features.rds" %>%
        sformat(list(name = filename))
      )
    )
    print("* Finished")
  } else{
    print("** Warning: Be caregul. Load hclust analysis")
    hclust_features <- readRDS(
      file.path(
        "..",
        "RDS",
        "{name}/dndrgm/single_linkage_features.rds" %>%
          sformat(list(name = filename))
      )
    )
    print("* Finished")
  }
  print(
    hclust_features
  )
}
# HEAD ----
folder <- "merged"
distances <- "tracto2016"
model <- "normal_GB_GB"
csv.name <- "fln_norm"
labels.name <- "imputation_labels"
common.path <- paste(distances, model, sep = "/")
csv.path <- paste(
  folder, "imputation", common.path, paste0(csv.name, ".csv"),
  sep = "/"
)
labels.path <- paste(
  folder, "labels", common.path, paste0(labels.name, ".csv"),
  sep = "/"
)
plot.path <- "../plots"
path.list <- list(
  csv = csv.path,
  labels = labels.path,
  plot = plot.path,
  common = common.path,
  folder = folder,
  model = model,
  distances = distances)
main(path.list)
