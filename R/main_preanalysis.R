make.aks.matrix <- function(
  AKS, bourder, colors, labels, regions, inst,
  preffix="", suffix="", on=T) {
  if (on) {
    print("** Plotting Aks matrix")
    source("functions/plot_aks_matrix.R")
    plot.aks.matrix(
      AKS, bourder, colors, labels, regions,
      path = inst$plot, foldername = inst$mainfolder,
      subfolder = "Preprocessing", preffix = preffix, suffix = suffix
    )
  } else
    print("** No Aks matrix")
}

make.aki.aik.density <- function(
  AKI, AIK, labels, regions, inst,
  suffix="", on=T) {
  if (on) {
    print("** Plotting Aks density")
    source("functions/plot_aki_aik_density.R")
    plot.aki.aik.density(
      AKI, AIK, labels, regions,
      path = inst$plot, foldername = inst$mainfolder,
      subfolder = "Preprocessing", suffix = suffix
    )
  } else
    print("** No Aks density")
}

make.aki.aik.scatter <- function(
  net, AKI, AIK, nt, labels, inst,
  suffix="", on=T) {
  if (on) {
    print("** Plotting scatter plots")
    source("functions/plot_aki_aik_scatter.R")
    plot.aki.aik.scatter(
      net, AKI, AIK, nt, labels,
      path = inst$plot, foldername = inst$mainfolder,
      subfolder = "Preprocessing", suffix = suffix
    )
  } else
    print("** No scatter plots")
}

make.fln.histogram <- function(net, inst, on=T) {
  if (on) {
    print("Plotting fln histogram")
    source("functions/plot_fln_histogram.R")
    plot.fln.histogram(
      net,
      path = inst$plot, foldername = inst$mainfolder,
      subfolder = "Preprocessing"
    )
  } else
    print("No fln histogram")
}

make.performance <- function(df, labels, inst, on=T) {
  if (on) {
    print("** Plotting performance")
    source("functions/plot_performance.R")
    plot.performance(
      df$net, df$aik.supra, df$aik.infra, df$aik, df$supra, df$infra, labels,
      path = inst$plot, foldername = inst$mainfolder,
      subfolder = "Preprocessing"
    )
  } else
    print("**No performance")
}

make.features <- function(
  net, distance, net.aik, net.aki, nt, nodes, names, inst, on=T
  ) {
  if (on) {
    source("functions/plot_features.R")
    print("** Plotting features")
    net[net == 0] <- NA
    distance <- distance %>%
      df.to.adj() %>%
      standardized.block.diag(nt)
    net.aik <- net.aik %>%
      standardized.diag(nt)
    net.aki <- net.aki %>%
      standardized.diag(nt)
    plot.features(
      nt, names,
      net, distance, net.aik, net.aki,
      path = inst$plot,
      folder = inst$folder,
      subfolder = "PreAnalysis"
    )
  } else
    print("** No features")
}

col.sum <- function(A, N) {
  x <- c()
  for ( i in 1:N) {
    x <- c(x, sum(A[, i], na.rm = T))
  }
  return(x)
}

get.mean.similarity <- function(df, ntrn, ntgt, nsrc=107) {
  net <- df$net
  AIK <- df$aik
  AKI <- df$aki
  AIK[diag(nsrc) == 1] <- NA
  AKI[diag(ntrn) == 1] <- NA
  similarity <- matrix(0, nrow = nsrc, ncol = ntgt)
  for (i in 1:ntrn) {
    for (j in 1:ntrn) {
      if (i != j) {
        jacp.aik <- AIK[i, j]
        jacp.aki <- AKI[i, j]
        similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
      }
    }
  }
  for (i in (ntrn + 1):nsrc) {
    for (j in 1:ntrn) {
      Ni <- which(net[i, ] > 0)
      Nj <- which(net[j, ] > 0)
      NiNj <- Ni[Ni %in% Nj]
      jacp.aik <- AIK[i, j]
      jacp.aki <- mean(AKI[j, NiNj], na.rm = T)
      similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
    }
  }
  return(similarity)
}

get.datasets <-  function(net, labels, nt, inst, filename="") {
  source("functions/adj_to_df.R")
  source("functions/df_to_adj.R")
  # Ensemble data ----
  ## Get data ----
  net <- net[net$target <= nt, ] %>%
    df.to.adj() %>%
    adj.to.df()
  source("functions/get_tracto2016.R")
  distance <- get.tracto2016(labels)
  distance <- distance[, 1:nt] %>%
    adj.to.df()
  return(
    list(
      net = net,
      distance = distance
    )
  )
}

transform.datasets <- function(dataset, nodes, nt, inst, filename="") {
  source("functions/compute_aik.R")
  source("functions/compute_aki.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/standardized.R")
  # Define datasets ----
  net <- dataset$net  %>%
    df.to.adj()
  distance <- dataset$distance
  # Compute Aks ----
  net.aik <- net %>%
    adj.to.df() %>%
    compute.aik(nt)
  net.aki <- net %>%
    adj.to.df() %>%
    compute.aki(nt)
  # Plot features ----
  make.features(
    net, distance, net.aik, net.aki, nt, nodes,
    c("w_standard", "dist", "jaccp_src", "jaccp_tgt"), inst,
    on = T
  )
}

compare_target_similarity <- function(dataset, nodes, nt, inst, filename="") {
  source("functions/compute_aki.R")
  source("functions/compute_aki_ec.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/standardized.R")
  # Define datasets ----
  net <- dataset$net  %>%
    df.to.adj()
  # Compute Aks ----
  net_aki <- net %>%
    adj.to.df() %>%
    compute.aki(nt) %>%
    adj.to.df() %>%
    dplyr::filter(source < nt)
  net_aki_ec <- net %>%
    adj.to.df() %>%
    compute.aki.ec(nt) %>%
    adj.to.df() %>%
    dplyr::filter(source < nt)
  # Compute error ----
  tb <- dplyr::tibble(
    error = net_aki$weight - net_aki_ec$weight
  )
  print(
    shapiro.test(tb$error)
  )
  # Compute statistics ----
  tb_stat <- tb %>%
    dplyr::summarise(
      sd = sd(error) %>%
      round(3),
      mean = mean(error) %>%
      round(3)
    )
  # Plot ----
  p <- tb  %>%
  ggplot2::ggplot(
    ggplot2::aes(error)
  ) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = ..density..),
    color = "gray",
    bins = 30
  ) +
  ggplot2::xlab("Error") +
  ggplot2::geom_text(
    ggplot2::aes(x = -0.1, y = 10),
    label = "mean = %.3f, sd = %.3f" %>%
      sprintf(
        tb_stat$mean,
        tb_stat$sd
      )
  ) +
  ggplot2::theme_classic()
  return(p)
}

target_similarity_source <- function(dataset, nodes, nt, inst, filename="") {
  source("functions/compute_aki.R")
  source("functions/compute_aik.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/standardized.R")
  # Define datasets ----
  net <- dataset$net  %>%
    df.to.adj()
  # Compute Aks ----
  net_target_sim <- net %>%
    adj.to.df() %>%
    compute.aki(nt) %>%
    adj.to.df() %>%
    dplyr::filter(source < nt)
  net_source_sim <- net %>%
    adj.to.df() %>%
    compute.aik(nt) %>%
    adj.to.df() %>%
    dplyr::filter(source < nt)
  # Compute error ----
  p <- dplyr::tibble(
    target = net_target_sim$weight,
    source = net_source_sim$weight
  )  %>%
  ggplot2::ggplot(
    ggplot2::aes(source, target)
  ) +
  ggplot2::geom_point(size = 0.5) +
  ggplot2::xlab("Source similarity") +
  ggplot2::ylab("Target similarity") +
  ggplot2::geom_smooth(
    method = "lm",
    se = F
  ) +
  ggpubr::stat_cor(
    ggplot2::aes(
      label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")
    ),
    label.x.npc = 0.5
  ) +
  ggplot2::theme_classic()
  return(p)
}

target_similarity_argument <- function(p, q, inst) {
  plot <- cowplot::plot_grid(
    p, q, ncol = 2, nrow = 1, labels = "AUTO"
  )
  png(
    "%s/%s/%s/argument_target_similarity.png" %>%
      sprintf(inst$plot, inst$folder, "PreAnalysis"),
    width = 12, height = 5, units = "in", res = 200
  )
  print(plot)
  dev.off()
}

main <- function(inst) {
  library(magrittr)
  source("functions/adj_to_df.R")
  nlog10 <- T
  nt <- 50
  nodes <- 107
  inst$mainfolder <- paste(inst$folder, "PreAnalysis", sep = "/")
  dir.create(
    sprintf("%s/%s", inst$plot, inst$mainfolder),
    showWarnings = F
  )
  # Load data ----
  source("functions/load_net.R")
  netx <- load.net(inst)
  net <- netx$net
  labels <- netx$labels
  source("functions/format_labels.R")
  labels <- labels %>%
    format.labels()
  if (nlog10) {
    net$weight[net$weight != 0] <- log10(net$weight[net$weight != 0]) + 7
  }
  source("functions/get_regions_colors.R")
  cr <- get.regions.colors(labels, nt)
  # Analyzing data ----
  datasets <- get.datasets(net, labels, nt, inst, filename = "default_data")
  # transform.datasets(datasets, nodes, nt, inst, filename = "transform_data")
  # target_similarity_argument(
  #   compare_target_similarity(datasets, nodes, nt, inst),
  #   target_similarity_source(datasets, nodes, nt, inst),
  #   inst
  # )
  source("functions/measure_similarity_dispersion_b.R")
  measure_similarity_dispersion_b(net, nodes, nt, labels, inst)
  # Plotting ----
}

setwd("/Users/jmarti53/Documents/Projects/LINKCOMMUNITIES/MAC/220126/R")
### HEAD ----
folder <- "merged"
distances <- "original"
model <- "normal"
csv_name <- "fln"
labels_name <- "rlabels"  # imputation_labels  rlabels
# merged/imputation/tracto2016/zz_model/   91x40
common_path <- paste(distances, model, sep = "/")
csv_path <- paste(
  folder, "imputation", common_path, paste0(csv_name, ".csv"),
  sep = "/"
)
labels_path <- paste(
  folder, "labels", common_path, paste0(labels_name, ".csv"),
  sep = "/"
)
plot_path <- "../plots"
path_list <- list(csv = csv_path,
                  labels = labels_path,
                  plot = plot_path,
                  common = common_path,
                  folder = folder,
                  model = model,
                  distances = distances)

main(path_list)
### Plot link by areas by distance with color on the weights
### Change mean to root squared sum
### In and out similarities colored by weight