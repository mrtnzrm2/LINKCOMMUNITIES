plot.region.commship <-  function(
  net, hcluster, k, labels, memberships, regions,
  path="", foldername="", subfolder="", filename=""
){
  if (!pracma::strcmp(subfolder, "")) {
    dir.create(
      "%s/%s/%s" %>%
        sprintf(path, foldername, subfolder),
      showWarnings = FALSE
    )
  }
  source("functions/assign_commship.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/format_regions.R")
  print("** Warning: Be careful selecting the right regions")
  print("** Using hclust as passed")
  source('functions/assign_commship.R')
  net <- assign.commship(net, hcluster, k, 0)
  # print("*** Warning: You are using a reference hclust to munkres the current hclust")
  # source('functions/assign_commship_reference.R')
  # net <- assign.commship.reference(net, hcluster, k)
  # source('functions/format_lincomm.R')
  # net <- format.lincomm(net, hcluster, k)
  net$commship <- net$commship %>%
    as.factor()
  net$source.label <- labels[net$source]
  net$target.label <- labels[net$target]
  commships <- net$commship %>%
    unique() %>%
    sort()
  nm <- commships %>%
    length()
  area.count <- data.frame(nodes=rep(1:length(labels), nm), regions=rep(regions$REGION[match(labels, regions$AREA)], nm), areas=rep(labels, nm), commship=rep(commships, each=length(labels)), size=0)
  for (i in 1:nrow(net)) {
    area.count$size[which(area.count$nodes == net$source[i] & area.count$commship == net$commship[i])] <-
      area.count$size[which(area.count$nodes == net$source[i] & area.count$commship == net$commship[i])]  + 1
    area.count$size[which(area.count$nodes == net$target[i] & area.count$commship == net$commship[i])] <-
      area.count$size[which(area.count$nodes == net$target[i] & area.count$commship == net$commship[i])]  + 1
  }
  area.count$areas <- factor(
    area.count$areas,
    levels = regions$AREA[order(regions$REGION)]
  )
  p <- ggplot2::ggplot(area.count, ggplot2::aes(areas, size, fill = commship)) +
    ggplot2::facet_wrap(~regions, scales = "free") +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("Areas") +
    ggplot2::ylab("Number of links") +
    ggplot2::scale_fill_brewer(name = "LC", palette = "Set2") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = 5
      )
    )
  p_leg <- p %>%
    ggpubr::get_legend()
  # p <- p +
  #   ggplot2::theme(
  #     legend.position = "none"
  #   )
  q_color <- regions %>%
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
  q_leg <- q_color %>%
    ggpubr::get_legend()
  leg <- cowplot::plot_grid(
    p_leg, q_leg, nrow = 2, ncol = 1
  )
  pq <- cowplot::plot_grid(
    p, leg, rel_widths = c(9, 1)
  )
   png(
    "%s/%s/%s/%s.png" %>%
      sprintf(path, foldername, subfolder, filename),
    width = 13, height = 5, units = "in", res = 200
  )
  print(p)
  dev.off()
}

plot.wdis.commship <- function(net, hcluster, k, path="", foldername="", subfolder="", filename=""){
  
  if (!pracma::strcmp(subfolder, '')) {
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  source("functions/assign_commship.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  net <- assign.commship(net, hcluster, k, 0)
  nat <- df.to.adj(net)
  nat <- adj.to.df(nat)
  print("** Warning: Be careful choosing the right distance matrix. Be sure it has the same order as the fln matrix")
  source("functions/get_tracto2016.R")
  distances <- get.tracto2016(labels)
  distances <- distances %>%
    adj.to.df()
  net$commship <- net$commship %>%
    as.factor()
  net$distances <- distances$weight[nat$weight != 0]
  p <- ggplot2::ggplot(net, ggplot2::aes(distances, w, color = commship)) +
    ggplot2::geom_point(size=0.3) +
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::theme_classic()
   png(
    "%s/%s/%s/%s.png" %>%
      sprintf(path, foldername, subfolder, filename),
    width = 13, height = 5, units = "in", res = 200
  )
  print(pq)
  dev.off()
}

plot.area.commship <- function(
  net, hcluster, k, labels, memberships, regions,
  path="", foldername="", subfolder="", filename=""
) {
  if (!pracma::strcmp(subfolder, "")) {
    dir.create(
      "%s/%s/%s" %>%
        sprintf(path, foldername, subfolder),
      showWarnings = FALSE
    )
  }
  source("functions/assign_commship.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  print("** Using hclust as passed")
  source('functions/assign_commship.R')
  net <- assign.commship(net, hcluster, k, 0)
  # print("*** Warning: You are using a reference hclust to munkres the current hclust")
  # source('functions/assign_commship_reference.R')
  # net <- assign.commship.reference(net, hcluster, k)
  # source('functions/format_lincomm.R')
  # net <- format.lincomm(net, hcluster, k)
  net$commship <- net$commship %>%
    as.factor()
  net$source.label <- labels[net$source]
  net$target.label <- labels[net$target]
  commships <- net$commship %>%
    unique() %>%
    sort()
  nm <- commships %>%
    length()
  area.count <- data.frame(
    nodes = rep(1:length(labels), nm),
    nodes.cluster = rep(memberships, nm),
    areas = rep(labels, nm),
    commship = rep(commships, each = length(labels)),
    size = 0
  )
  for (i in 1:nrow(net)) {
    area.count$size[
      which(
        area.count$nodes == net$source[i] &
          area.count$commship == net$commship[i]
      )
    ] <- area.count$size[
      which(
        area.count$nodes == net$source[i] &
          area.count$commship == net$commship[i]
        )
      ] + 1
    area.count$size[
      which(
        area.count$nodes == net$target[i] &
          area.count$commship == net$commship[i]
      )
    ] <- area.count$size[
      which(
        area.count$nodes == net$target[i] &
          area.count$commship == net$commship[i]
      )
    ] + 1
  }
  tr.labels <- labels[order(memberships)]
  area.count$areas <- factor(area.count$areas, levels = tr.labels)
  p <- area.count %>%
    ggplot2::ggplot(
      ggplot2::aes(
        areas, size, fill = commship
      )
    ) +
    ggplot2::facet_wrap(~nodes.cluster, scales = "free") +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("Areas") +
    ggplot2::ylab("Number of links") +
    ggplot2::scale_fill_brewer(name = "LC", palette = "Set2") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        color = regions$COLOR[
          match(tr.labels, regions$AREA)
        ],
        angle = 90,
        size = 5
      )
    )
  p_leg <- p %>%
    ggpubr::get_legend()
  p <- p +
    ggplot2::theme(
      legend.position = "none"
    )
  q_color <- regions %>%
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
  q_leg <- q_color %>%
    ggpubr::get_legend()
  leg <- cowplot::plot_grid(
    p_leg, q_leg, nrow = 2, ncol = 1
  )
  pq <- cowplot::plot_grid(
    p, leg, rel_widths = c(9, 1)
  )
  png(
    "%s/%s/%s/%s.png" %>%
      sprintf(path, foldername, subfolder, filename),
    width = 13, height = 5, units = "in", res = 200
  )
  print(pq)
  dev.off()
}