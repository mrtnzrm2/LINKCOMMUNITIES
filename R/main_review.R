load.merde <- function(
  process, nodes, labels, net, net_cluster,
  leaves, filename, save=T
  ) {
  if (save) {
    print("Loading the merde")
    dir.create(sprintf("../RDS/%s/merde/", filename), showWarnings = F)
    source("functions/la_arbre_a_merde.R")
    ham_merde <- la.arbre.a.merde(
      process, nodes, labels[1:nodes],
      net, net_cluster, leaves
    )
    saveRDS(
      ham_merde,
      paste(
        "../RDS", filename, "merde", "node_agglomeration.rds", sep = "/"
      )
    )
  } else{
    print("Not loading the merde")
    ham_merde <- readRDS(
      paste(
        "../RDS", filename, "merde", "node_agglomeration.rds", sep = "/"
      )
    )
  }
  return(ham_merde)
}

make.cascade <- function(net, merde, nodes, labels, regions, inst, on=T) {
  if (on) {
    print("Plotting the merde cascade")
    source("functions/cascade_sil_vous_plaît.R")
    cascade_merde <- cascade.sil.vous.plaît(net, merde, nodes)
    source("functions/cascade_arc_en_ciel.R")
    cascade.arc.en.ciel(
      merde, nodes, labels[1:nodes], regions, cascade_merde,
      path = inst$plot,
      subfolder = "ArcInCiel",
      plt = T,
      filename = "Arco_byH",
      foldername = inst$mainfolder
    )
  } else
    print("No cascade")
}

make.heuristics <- function(features, inst, on=T) {
  if (on) {
    print("Printing parameters plots")
    source("functions/plot_heuristics.R")
    plot.heuristics(
      inst$plot,
      inst$mainfolder,
      features,
      subfolder = "Heuristics"
    )
  } else
    print("No parameters")
}

make.mountain <- function(merde, labels, nodes, coords, regions, inst, on=T) {
  if (on) {
    print("Plotting the mountain view")
    source("functions/mountain_view.R")
    mountain.view(
      merde,
      labels = labels[1:nodes],
      net.coord = coords,
      regions = regions,
      path = inst$plot,
      nodes = nodes,
      subfolder = "Mountain",
      filename = "ContourPlot",
      foldername = inst$mainfolder,
      plt = T
    )
  } else
    print("No mountain view")
}

make.circular <- function(Ns, features, merde, inst, regions, on=T){
  if (on) {
    print("Plotting circular dendrogram")
    source("functions/plot_dendrogram_circular.R")
    for (k.m.2 in Ns) {
      k.2 <- features$K[which(features$NEC == k.m.2)]
      if (length(k.2) == 0)
        next
      k.2 <- k.2[length(k.2)]
      plot.dendrogram.circular(
        merde,
        k = k.m.2,
        filename = sprintf("polarMerde_k_%i", k.m.2),
        plt = T,
        alternative = F,
        tip.labels.size = 3.5,
        line.size = 1,
        legend.plot.size = 30, strip.text.size = 1,
        offset.strip.text = 0.05, offset.strip = 0.5,
        branch.none = T,
        polar.text.size = 5, bar.line.size = 3,
        subfolder = "BrainEvolutionSim",
        path = inst$plot,
        tree.regions = list(T, regions),
        foldername = inst$mainfolder
      )
    }
  } else
    print("No circular dendrogram")
}

make.dend.matrix <- function(
  Ns, features, merde, cluster, nodes, inst, regions, on=T
  ) {
  if (on) {
    print("Plotting matrix dendrogram")
    source("functions/plot_dendrogram_matrix.R")
    for (k.m.2 in Ns) {
      k.2 <- features$K[which(features$NAC == k.m.2)]
      if (length(k.2) == 0)
        next
      k.2 <- k.2[length(k.2)]
      plot.dendrogram.matrix(
        net[net$source <= nodes & net$target <= nodes, ],
        merde, cluster, labels[1:nodes], regions, nodes,
        k = k.2,
        height = 0,
        k.merde = k.m.2,
        filename = sprintf("recMerde_k_%i_v2", k.m.2),
        foldername = inst$mainfolder,
        plt = T,
        alternative = F,
        line.size = 1,
        tip.labels.size = 1,
        branch.none = F,
        subfolder = "Matrices_ns",
        path = inst$plot
      )
    }
  } else
    print("No matrix dendrogram")
}

make.link.comm <- function(
  Ns, net, nodes, labels, cluster, leaves,
  regions, memberships, inst, on=T
  ) {
  if (on) {
    print("Plotting link communities matrices")
    source("functions/plot_linkcomm_matrix.R")
    for (k.2 in Ns) {
      plot.lincomm.matrix(
        net[net$source <= nodes & net$target <= nodes, ],
        cluster,leaves, nodes, k.2, 0,
        labels, regions,
        memberships = memberships,
        plt = T,
        subfolder = "Matrices",
        path = inst$plot,
        filename = "%s_linkcommunities_%s" %>%
          sprintf(k.2, inst$suffix),
        foldername = inst$mainfolder)
    }
  } else
    print("No link communities matrix")
}

make.link.comm_b <- function(
  Ns, net, nodes, labels, cluster, leaves,
  regions, memberships, inst, on=T
  ) {
  if (on) {
    print("Plotting link communities matrices b-style")
    source("functions/plot_linkcomm_matrix_b.R")
    for (k.2 in Ns) {
      plot.lincomm.matrix_b(
        net[net$source <= nodes & net$target <= nodes, ],
        cluster,leaves, nodes, k.2, 0,
        labels, regions,
        memberships = memberships,
        plt = T,
        subfolder = "Matrices",
        path = inst$plot,
        filename = "%s_linkcommunities_%s" %>%
          sprintf(k.2, inst$suffix),
        foldername = inst$mainfolder)
    }
  } else
    print("No link communities matrix")
}

make.heatmap.communities <- function(
  Ns, net, nodes, labels, cluster, leaves,
  regions, memberships, inst, on=T
  ) {
  if (on) {
    print("Plotting link communities matrices")
    source("functions/plot_heatmap_communities.R")
    for (k.2 in Ns) {
      plot.heatmap.communities(
        net[net$source <= nodes & net$target <= nodes, ],
        cluster,leaves, nodes, k.2, 0,
        labels, regions,
        memberships = memberships,
        plt = T,
        subfolder = "HeatMap",
        path = inst$plot,
        filename = "%s_heatmap_%s" %>%
          sprintf(k.2, inst$suffix),
        foldername = inst$mainfolder)
    }
  } else
    print("No link communities matrix")
}

make.heatmap.communities_b <- function(
  Ns, net, nodes, labels, cluster, leaves,
  regions, memberships, inst, on=T
  ) {
  if (on) {
    print("Plotting link communities matrices with b-style")
    source("functions/plot_heatmap_communities_b.R")
    for (k.2 in Ns) {
      plot.heatmap.communities_b(
        net[net$source <= nodes & net$target <= nodes, ],
        cluster,leaves, nodes, k.2, 0,
        labels, regions,
        memberships = memberships,
        plt = T,
        subfolder = "HeatMap",
        path = inst$plot,
        filename = "%s_heatmap_%s" %>%
          sprintf(k.2, inst$suffix),
        foldername = inst$mainfolder)
    }
  } else
    print("No link communities matrix")
}

make.bars <-  function(
  Ks, net, merde, cluster, labels,
  regions, inst, on=T) {
  if (on) {
    print("Plotting link comms bar plots")
    source("functions/plot_bar_link_comms.R")
    for (k in Ks)
      plot.bar.link.comms(
        k, net, merde, cluster, labels, regions,
        path = inst$plot,
        foldername = inst$mainfolder,
        subfolder = "BarLinkComms",
        filename = paste0(k, "_bar")
      )
  } else
    print("No bars")
}

save.tables <- function(Ks, net, merde, labels, inst, on=T) {
  if (on) {
    print("Creating tables")
    source("functions/link_table.R")
    for (k in Ks) {
      dir.create(
        paste(
          "../CSV", inst$folder, "tables", inst$common, k, sep = "/"
        ),
        showWarnings = F
      )
      memberships <- link.table(net, merde, labels, k)
      h1 <- memberships[memberships$DIC == "H", ]
      v1 <- memberships[memberships$DIC == "V", ]
      commships <- memberships$COMMSHIP %>%
        unique() %>%
        sort()
      nc <- commships %>%
        length()
      rc <- 1:nc
      for (j in 1:nc) {
        h2 <- h1[h1$COMMSHIP == commships[j], ]
        if (nrow(h2) < length(labels)) {
          h2 <- h2 %>%
            rbind(
              data.frame(
                AREA = labels[!labels %in% h2$AREA],
                COMMSHIP = commships[j],
                P = 0,
                DIC = "H"
              )
            )
        }
        h2 <- h2[match(labels, h2$AREA), ]
        v2 <- v1[v1$COMMSHIP == commships[j], ]
        if (nrow(v2) < length(labels)) {
          v2 <- v2 %>%
            rbind(
              data.frame(
                AREA = labels[!labels %in% v2$AREA],
                COMMSHIP = commships[j],
                P = 0,
                DIC = "V"
              )
            )
        }
        v2 <- v2[match(labels, v2$AREA), ]
        write.csv(
          h2,
          paste(
            "../CSV", inst$folder, "tables",
            inst$common, k, paste0("h_", rc[j], ".csv"),
            sep = "/"
          ),
          row.names = F
        )
        write.csv(
          v2,
          paste(
            "../CSV", inst$folder, "tables",
            inst$common, k, paste0("v_", rc[j], ".csv"),
            sep = "/"
          ),
          row.names = F
        )
      }
      max.membership <- memberships %>%
        dplyr::group_by(AREA) %>%
        dplyr::summarise(
          COMMSHIP = COMMSHIP[which(P == max(P))],
          P = max(P)
        )
      max.membership <- max.membership[match(labels, max.membership$AREA), ]
      write.csv(
        max.membership,
        paste(
          "../CSV", inst$folder, "tables", inst$common, k, "max.csv", sep = "/"
        ),
        row.names = F
      )
    }
  } else
    print("No creating tables")
}

save.tables.wsbm <- function(K, labels, memberships, inst, extra="", on=T) {
  if (on) {
    print("Saving WSBM table")
    print("----> Warning: Choose the right partition")
    for (k in K) {
      "../CSV/%s/tables/%s/%i" %>%
        sprintf(inst$folder, inst$common, k) %>%
        dir.create(showWarnings = F, recursive = T)
      memberships <- memberships %>%
        unlist() %>%
        as.numeric()
      dplyr::tibble(
        AREA = labels,
        COMMSHIP = memberships,
        DIR = "WSBM") %>%
        write.csv(
          "../CSV/%s/tables/%s/%i/%s%s.csv" %>%
          sprintf(
            inst$folder, inst$common, k, inst$suffix, extra
          )
        )
    }
  } else{
    print("No WSBM table")
  }
}

lc_color_area <- function(labels, memberships) {
  color_lc_areas <- rep(0, length(labels))
  color_lc_areas[
    memberships == 1 |
    memberships == 2
  ] <- "#fc8d62"
  color_lc_areas[
    memberships == 4 |
    memberships == 5
  ] <- "#8da0cb"
  # color_lc_areas[
  #   memberships == 3
  # ] <- "#e78ac3"
  color_lc_areas[
    memberships == 3
  ] <- "#fdbf6f"
  color_lc_areas[
    memberships == 6
  ] <- "#66c2a5"
  return(color_lc_areas)
}

save.clc <- function(
  K, labels, memberships,
  inst, on=T) {
  if (on) {
    print("Saving WSBM table")
    print("----> Warning: Choose the right partition")
    for (k in K) {
      memberships <- memberships %>%
        unlist() %>%
        as.numeric()
      dplyr::tibble(
        AREA = labels,
        COMMSHIP = memberships,
        LC = lc_color_area(labels, memberships),
        P = 0.7,
        DIR = "WSBM") %>%
        write.csv(
          "../CSV/%s/tables/%s/%i/clc_%s.csv" %>%
            sprintf(inst$folder, inst$common, k, inst$suffix)
        )
    }
  } else{
    print("No WSBM table")
  }
}

make.hier <- function(
  K, net, net.cluster, regions, labels, membership,
  inst, on=T) {
  if (on) {
    print("Printing Hierarchical Edge Bundle")
    source("functions/plot_hier_edg_bundle.R")
    for (k in K) {
      plot.hier.edg.bundle(
        k, net, net.cluster, regions, labels, membership,
        path = inst$plot,
        foldername = inst$mainfolder,
        subfolder = "HEB",
        filename = paste(k, inst$suffix)
      )
    }
  } else{
    print("No hierarchical edge bundle")
  }
}

make.bar.block.sim <- function(R, K, net, net.cluster, inst, on=T) {
  if (on) {
    print("Plotting bar-block-sim-v1")
    source("functions/plot_bar_block_sim.R")
    for (r in R) {
      for (k in K) {
        plot.bar.block.sim(
          r, k, net, net.cluster,
          path = inst$plot,
          foldername = inst$mainfolder,
          subfolder = "BlockSim",
          filename = sprintf("k_%i_r_%i", k, r)
        )
      }
    }
  } else
    print("No bar-block-sim_v1")

}

make.densities <- function(net, labels, membership, on=T) {
  if (on) {
    print("Plotting block densities")
    source("functions/plot_densities.R")
    plot.densities(net, labels, membership)
  } else
    print("No densities")
}

make.Astats.commship <- function(net, hcluster, k, labels, memberships, regions,
inst, on=T) {
  if (on) {
    source("functions/plot_stats_commship.R")
    print("* Plotting different statistical plots")
    plot.area.commship(
      net, hcluster, k, labels, memberships, regions,
      path = inst$plot,
      foldername = inst$mainfolder,
      subfolder = "stats",
      filename = paste("area_size", "commship", k, inst$suffix, sep = "_")
    )
  } else
    print("No Astat commship")
}

make.Rstats.commship <- function(net, hcluster, k, labels, memberships, regions,
inst, on=T) {
  if (on) {
    source("functions/plot_stats_commship.R")
    print("* Creating different statistical plots")
    plot.region.commship(
      net, hcluster, k, labels, memberships, regions,
      path = inst$plot,
      foldername = inst$mainfolder,
      subfolder = "stats",
      filename = paste("region_size", "commship", k, inst$suffix, sep = "_")
    )
  } else
    print("No Rstat commship")
}

save.lincomm.matrix <- function(net, hcluster, k, inst, save=T) {
  if (save) {
    print("* Saving Commship matrix in WSBM folder")
    dir.create(
      paste(
        "../WSBM/CD/CSV/commships",
        inst$folder,
        inst$common,
        sep = "/"
      ),
      showWarnings = F
    )
    # source("functions/format_lincomm.R")
    # net <- format.lincomm(net, hcluster, k)
    source("functions/assign_commship.R")
    net <- assign.commship(net, hcluster, k, 0)
    source("functions/df_to_adj.R")
    net <- with(
      net,
      data.frame(
        source = source,
        target = target,
        weight = commship
      )
    )
    net <- net %>%
      df.to.adj()
    write.csv(
      net,
      paste(
        "../WSBM/CD/CSV/commships",
        inst$folder,
        inst$common,
        "%i_commships.csv" %>%
          sprintf(k),
        sep = "/"
      ),
      row.names = F
    )
  } else
    print("Not saving link community matrix")
}

make.area.hierarchy <- function(
  area, K, leaves, nodes, labels,
  net, filename, inst, save=T
) {
  if (save) {
    if (is.character(area)) {
      area <- which(labels == area)
      print(
        "* Start computing area history: %s" %>%
          sprintf(labels[area])
      )
    } else if (is.numeric(area)) {
      print(
        "* Start computing area history: %s" %>%
          sprintf(labels[area])
      )
    }
    source("functions/adj_to_df.R")
    source("functions/df_to_adj.R")
    in_deg <- net %>%
      df.to.adj()
    in_deg[in_deg > 0] <- 1
    out_deg <- in_deg %>%
      apply(1, sum)
    in_deg <- in_deg %>%
      apply(2, sum)
     net_sim <- readRDS(
      "../RDS/%s/similarity/link_similarity.rds" %>%
        sprintf(filename)
    )
    net_sim[net_sim == -1] <- NA
    net_dis <- -log10(net_sim)
    net_dis[is.na(net_dis)] <- max(net_dis, na.rm = T) + 1
    Rcpp::sourceCpp("../cpp/area_hierarchy.cpp")
    source("functions/fx.R")
    hierarchy <- area_hierarchy(
      area, leaves, nodes,
      net$source, net$target,
      in_deg, out_deg,
      labels,
      fx(net_dis)
    )
    labels_hierarchy <- labels[hierarchy[hierarchy[, 1] > 0, 1]]
    K_hierarchy <- hierarchy[hierarchy[, 1] > 0, 2] %>%
      log10()
    K_hierarchy <- K_hierarchy / max(K_hierarchy)
    K_hierarchy <- c(
      K_hierarchy,
      rep(-1, length(labels) - length(labels_hierarchy))
    )
    labels_hierarchy <- c(
      labels_hierarchy,
      labels[!labels %in% labels_hierarchy]
    )
    # Has to have same order as labels
    labels_order <- labels %>%
      match(labels_hierarchy)
    labels_hierarchy <- labels_hierarchy[labels_order]
    K_hierarchy <- K_hierarchy[labels_order]
    save.tables.wsbm(
      K, labels_hierarchy, K_hierarchy, inst,
      extra = "_local_%s" %>% sprintf(labels[area]), on = T
    )
  } else {
    print("* Not area hierarchy")
    print(which(labels == labels[area]))
    print(labels[which(labels == labels[area])])
  }
}

node_commships_format <- function(inst) {
  ### Load WSBM partition
  print("Warning: Be careful choosing the right partition")
  node_memberships <- read.csv(inst$wsbm, header = F) %>%
    unlist() %>%
    as.numeric()
  # Map the membership to the interval 1-K
  # This is good if there are missing comms in
  # the file.
  source("functions/arrange_memberships_values.R")
  node_memberships <- node_memberships %>%
    arrange.memberships.values()
  # # Compute munkers between the selected labels
  # # and one that is for reference inside the function.
  # print("** Arrange memberships with respect to a reference")
  # source("functions/arrange_memberships_reference.R")
  # node_memberships <- arrange.memberships.reference(
  #   node_memberships,
  #   inst$R,
  #   inst
  # )
  # Once nodes' arrengement is clear, you can
  # permute some of them to make it look better.
  print("** Arrange membership with some manual choices")
  source("functions/arrange_membership_manual.R")
  # node_memberships <- node_memberships %>%
  #   arrange.membership.manual(
  #     c(
  #       "1" = 5,
  #       "3" = 6
  #     )
  #   )
  ###
  node_memberships <- node_memberships %>%
    arrange.membership.manual(
      c(
        "2" = 5,
        "6" = 3
      )
    )
  node_memberships <- node_memberships %>%
    arrange.membership.manual(
      c(
        "5" = 7,
        "1" = 2
      )
    )
  node_memberships <- node_memberships %>%
    arrange.membership.manual(
      c(
        "2" = 3,
        "5" = 6
      )
    )
  return(node_memberships)
}

main <- function(inst) {
  linkage <- "average"
  index <- "jaccp"
  mode <- "ALPHA"
  nlog10 <- T
  EC <- F
  tag <- paste(inst$model, inst$distances, inst$folder, sep = "_")
  source("functions/model_name_analysis.R")
  filename <- model.name(index, linkage, mode, nlog10, EC, tag)
  source("functions/sformat.R")
  inst$mainfolder <- paste(
    inst$folder,
    inst$common,
    paste(toupper(linkage), "full", "l10", sep = "_"),
    sep = "/"
  )
  dir.create(sprintf("%s/%s", inst$plot, inst$mainfolder), showWarnings = F)
  #### Load network
  source("functions/load_net.R")
  netx <- load.net(inst)
  net <- netx$net
  nodes <- netx$nodes
  leaves <- netx$leaves
  labels <- netx$labels
  coords <- netx$fln.coords
  regions <- netx$regions
  if (nlog10) {
    net$weight[net$weight != 0] <- -log10(net$weight[net$weight != 0])
  }
  # Node memberships ----
  node_memberships <- node_commships_format(inst)
  ### Load data
  net.cluster <- readRDS(
    "../RDS/%s/hclust/hierarchical_clustering.rds" %>%
    sprintf(filename)
  )
  hclust.features <- readRDS(
    "../RDS/%s/dndrgm/single_linkage_features.rds" %>%
    sprintf(filename)
  )
  # ham.merde <- load.merde(
  #   hclust.features,
  #   nodes,
  #   labels,
  #   net,
  #   net.cluster,
  #   leaves,
  #   filename, save = T
  # )
  ### Plots start here:
  make.heuristics(hclust.features, inst, on = F)
  make.area.hierarchy(
    "v1pcuf", 4, leaves, nodes, labels,
    net, filename, inst, save = T
  )
  make.cascade(net, ham.merde, nodes, labels, regions, inst, on = F)
  make.mountain(ham.merde, labels, nodes, coords, regions, inst, on = F)
  make.circular(4, hclust.features, ham.merde, inst, regions, on = F)
  make.dend.matrix(
    Ns, features, ham.merde, net.cluster, nodes,
    inst, regions, on = F
  )
  make.link.comm(
    2:7, net, nodes, labels, net.cluster, leaves,
    regions, node_memberships, inst, on = F
  )
  make.link.comm_b(
    4, net, nodes, labels, net.cluster, leaves,
    regions, node_memberships, inst, on = F
  )
  make.heatmap.communities(
    4, net, nodes, labels, net.cluster, leaves,
    regions, node_memberships, inst, on = F
  )
  make.heatmap.communities_b(
    4, net, nodes, labels, net.cluster, leaves,
    regions, node_memberships, inst, on = F
  )
  make.bars(107:2, net, ham.merde, net.cluster, labels, regions, inst, on = F)
  make.hier(
    4, net, net.cluster, regions, labels, node_memberships, inst,
    on = F
  )
  make.bar.block.sim(6, 4, net, net.cluster, inst, on = F)
  make.densities(net, labels, node_memberships, on = F)
  make.Astats.commship(
    net, net.cluster, 2, labels, node_memberships, regions,
    inst, on = F
  )
  make.Rstats.commship(
    net, net.cluster, 2, labels, node_memberships, regions,
    inst, on = F
  )
  ## Table part:
  save.tables(34, net, net.cluster, labels, inst, on = F)
  save.tables.wsbm(4, labels, node_memberships, inst, on = F)
  save.lincomm.matrix(net, net.cluster, 3, inst, save = F)
  save.clc(3, labels, node_memberships, inst, on = F)

}

### HEAD ###
library(magrittr)
data <- "commships"
folder <- "merged"
distances <- "tracto2016"
model <- "normal_GB_GB"
csv.name <- "fln_norm"
MAN <- "MAX"
R <- 7
al <- 0
commships <- 4
labels.name <- "imputation_labels"
wsbm.name <- "al_%i/%i_LCS/K/%s/k_%i.csv" %>%
  sprintf(al, commships, MAN, R)
suffix <- "wsbm_al_%i_%s_k_%i" %>%
  sprintf(al, MAN, R)
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
wsbm.path <- paste(
  "../WSBM/CD/CSV/labels", data, folder, common.path, wsbm.name,
  sep = "/"
)

path.list <- list(
  csv = csv.path,
  labels = labels.path,
  plot = plot.path,
  common = common.path,
  folder = folder,
  model = model,
  distances = distances,
  wsbm = wsbm.path,
  R = R,
  al = al,
  MAN = MAN,
  suffix = suffix
)

main(path.list)