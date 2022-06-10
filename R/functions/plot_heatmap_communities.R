label_axis_position_x <- function(memberships) {
  get.order <- order(memberships)
  memberships <- memberships[get.order]
  nm <- memberships %>%
    unique() %>%
    length()
  tab.mem <- memberships %>%
    table()
  bourders <- rep("", nm - 1)
  for (i in 1:nm) {
    bourders[i] <- sum(tab.mem[1:i])
  }
  bourders <- bourders[1:(length(bourders) - 1)] %>%
    as.numeric()
  return(bourders)
}

label_axis_position_y <- function(memberships, label) {
  get.order <- order(memberships)
  memberships <- memberships[get.order]
  nm <- memberships %>%
    unique() %>%
    length()
  tab.mem <- memberships %>%
    table()
  bourders <- rep("", nm - 1)
  for (i in 1:nm) {
    bourders[i] <- sum(tab.mem[1:i])
  }
  bourders <- bourders[1:(length(bourders) - 1)] %>%
    as.numeric()
  bourders <- length(label) - bourders
  return(bourders)
}

plot.heatmap.communities <- function(
  net, net.cluster,
  leaves, nodes, k, best.height, labels,
  regions, memberships, plt=F,
  filename="", foldername="", subfolder="",
  path="") {
  source("functions/extract_k_partition.R")
  source("functions/with_potholes.R")
  source("functions/linkcomm_parameters.R")
  source("functions/gg_color_hue.R")
  source("functions/get_bourder_labels.R")
  source("functions/arrange_membership_manual.R")
  if (!pracma::strcmp(subfolder, "")) {
    dir.create(
      "%s/%s/%s" %>%
      sprintf(
        path, foldername, subfolder
      ),
      showWarnings = FALSE
    )
  }
  print("** Using hclust as passed")
  source('functions/assign_commship.R')
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  net <- assign.commship(net, net.cluster, k, 0)
  # print("*** Warning: You are using a reference hclust to munkres the current hclust")
  # source("functions/format_lincomm.R")
  # net <- format.lincomm(net, net.cluster, k)
  bourders_x <- label_axis_position_x(memberships)
  bourders_y <- label_axis_position_y(memberships, labels)
  tr.labels <- labels[order(memberships)]
  A <- net %>%
    with(
      dplyr::tibble(
        source = source,
        target = target,
        weight = w
      )
    ) %>%
    df.to.adj() %>%
    adj.to.df()
  A$weight[A$weight == 0] <- NA
  A$source_label <- labels[A$source]
  A$target_label <- labels[A$target]
  A$source_label <- factor(A$source_label, levels = rev(tr.labels))
  A$target_label <- factor(A$target_label, levels = tr.labels)
  A$LC <- A$weight %>%
    as.factor()
  q <- A %>%
    ggplot2::ggplot(
      ggplot2::aes(
        target_label, source_label, fill = weight
      )
    ) +
    ggplot2::geom_tile() +
    ggplot2::geom_vline(
      xintercept = bourders_x + 0.5,
      color = "green",
      size = 1
    ) +
    ggplot2::geom_hline(
      yintercept = bourders_y + 0.5,
      color = "green",
      size = 1
    ) +
    ggplot2::labs(
      x = "Target areas",
      y = "Source Areas"
    ) +
    viridis::scale_fill_viridis(
      na.value = "gray",
      direction = 1,
      option = "A",
      name = "log10(FLN)+7"
    ) +
    ggplot2::theme_classic() +
        ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        color = regions$COLOR[match(tr.labels, regions$AREA)]
      ),
      axis.text.y = ggplot2::element_text(
        color = regions$COLOR[match(rev(tr.labels),
        regions$AREA)]
      )
    )
    q_leg <- q %>%
      ggpubr::get_legend()
    q <- q +
      ggplot2::theme(legend.position = "none")
    p_color <- regions %>%
      ggplot2::ggplot(
        ggplot2::aes(
          AREA, REGION, fill = REGION
        )
      ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(
        name = "Regions",
        values = unique(regions$COLOR[order(regions$REGION)])
      )
    p_leg <- p_color %>%
      ggpubr::get_legend()
    leg <- cowplot::plot_grid(
      q_leg, p_leg, nrow = 2, ncol = 1
    )
    pq <- cowplot::plot_grid(
      q, leg, rel_widths = c(9, 1)
    )
  if (plt) {
    png(
      "%s/%s/%s/%s.png" %>%
      sprintf(
        path, foldername, subfolder, filename
      ),
      width = 14, height = 11, units = "in", res = 200
    )
    print(pq)
    dev.off()
  }
  else{
    print(pq)
  }
}