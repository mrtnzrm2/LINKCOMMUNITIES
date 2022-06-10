make.WSBM.matrix <- function(net, labels, nodes, path_al, inst, on=T, plt=T) {
  Data <- data.frame()
  if (on) {
    print("* Scanning all k.csv in the MAX and CON folders for all alphas")
    source("functions/plot_WSBM_matrix.R")
    source("functions/with_potatoes.R")
    if (dir.exists(path_al)) {
      al.dirs <- unlist(list.dirs(path_al, recursive = F))
      for (al in al.dirs) {
        al. <- unlist(stringr::str_split(al, "al_"))[2]
        alpha <- as.numeric(al.)
        if (is.na(alpha))
          next
        print("** Need to create some folders in the plot folder")
        dir.create(
          sprintf("%s/%s/ALPHA", inst$plot, inst$mainfolder),
          showWarnings = F
        )
        dir.create(
          sprintf("%s/%s/ALPHA/al_%i", inst$plot, inst$mainfolder, alpha),
          showWarnings = F
        )
        sub.dirs <- unlist(list.files(file.path(al, "K")))
        for (sub in sub.dirs) {
          print("*** Need to create some subfolders")
          dir.create(
            "%s/%s/ALPHA/al_%i/%s" %>%
              sprintf(inst$plot, inst$mainfolder, alpha, sub),
            showWarnings = F
          )
          k.dirs <- unlist(list.files(file.path(al, "K", sub)))
          K <- matrix(0, ncol = (nodes + 1), nrow = length(k.dirs)) %>%
            as.data.frame()
          colnames(K) <- c("K", paste0("id.", seq_len(nodes)))
          y <- 1
          for (k.f in k.dirs) {
            k. <- unlist(stringr::str_split(k.f, "k_"))[2]
            k. <- unlist(stringr::str_split(k., ".csv"))[1]
            k <- as.numeric(k.)
            labels.comms <- read.csv(
              "%s/al_%i/K/%s/k_%i.csv" %>%
                sprintf(path_al, alpha, sub, k), header = F
              ) %>%
                unname() %>%
                as.matrix() %>%
                pracma::Reshape(nodes, 1)
            K[y, 1] <- k %>%
              as.numeric()
            K[y, 2:ncol(K)] <- labels.comms[, 1] %>%
              as.numeric()
            y <- y + 1
            size.comms <- with.potatoes(labels.comms)
            size.comms <- dplyr::tibble(
              size = size.comms[[2]],
              label = size.comms[[1]]
            )
            size.comms <- size.comms[
              order(size.comms$size, decreasing = T),
            ]
            size.comms$relabel <- seq_len(nrow(size.comms))
            for (e in seq_len(labels.comms))
              labels.comms[e] <- size.comms$relabel[
                which(size.comms$label == labels.comms[e])
              ]
            "**** Plotting the WSBM matrix for alpha: %i and k: %i" %>%
              sprintf(alpha, k) %>%
              print()
            if (plt)
              plot.WSBM.matrix(
                net, labels[1:nodes], labels.comms,
                filename = sprintf("k_%i", k), subfolder = "MatrixPatition",
                path = "%s/%s/ALPHA/al_%i/%s" %>%
                  sprintf(inst$plot, inst$mainfolder, alpha, sub)
              )
            Data <- Data %>%
              dplyr::bind_rows(
                dplyr::tibble(
                  alpha = alpha / 100,
                  k = k,
                  membership = labels.comms[, 1],
                  labels = labels[1:nodes],
                  # x=fln.coods[,1],
                  # y=fln.coods[,2],
                  type = sub
                )
              )
          }
        }
      }
    }
  } else
    print("No WSBM matrix")
  return(Data)
}
make.WSBM.stats <- function(path_al, inst, on=T) {
  if (on) {
    print("* Scanning and comparing different characteristics of data")
    print("** Creating folders")
    "%s/%s/ALPHA/SUMMARY" %>%
      sprintf(inst$plot, inst$mainfolder) %>%
      dir.create(showWarnings = F, recursive = T)
    path.plot <- sprintf("%s/%s/ALPHA/SUMMARY", inst$plot, inst$mainfolder)
    LOGEVCON.DF <- dplyr::tibble()
    LOGEVMAX.DF <- dplyr::tibble()
    NMI.DF <- dplyr::tibble()
    PARTITIONS.DF <- dplyr::tibble()
    print("** Warning: Pay attention to the K and alpha ranges")
    K <- seq(2,15)
    alpha <- c(0)
    for (al in alpha) {
      for (k in K) {
        if (
          file.exists(
            "%s/al_%i/LOGEV/k_%i/logev_cons.csv" %>%
            sprintf(path_al, al, k)
          )
        ) {
          logev.con <- read.csv(
            "%s/al_%i/LOGEV/k_%i/logev_cons.csv" %>%
              sprintf(path_al, al, k), header = F
            ) %>%
              as.numeric()
          LOGEVCON.DF <- LOGEVCON.DF %>%
            dplyr::bind_rows(
              dplyr::tibble(
                al = al / 100,
                k = k,
                x = logev.con,
                type = "CON"
              )
            )
        }
        if (
          file.exists(
            "%s/al_%i/LOGEV/k_%i/logev_max.csv" %>%
              sprintf( path_al, al, k)
            )
          ) {
          logev.max <- read.csv(
            "%s/al_%i/LOGEV/k_%i/logev_max.csv" %>%
              sprintf(path_al, al, k), header = F
          ) %>%
            as.numeric()
          LOGEVMAX.DF <- LOGEVMAX.DF %>%
            dplyr::bind_rows(
              dplyr::tibble(
                al = al / 100,
                k = k,
                x = logev.max,
                type = "MAX"
              )
            )
        }
        if (
          file.exists(
            "%s/al_%i/NMI/k_%i/nmi_max_cons.csv" %>%
              sprintf(path_al, al, k)
          )
        ) {
          nmi.cos.max <- read.csv(
            "%s/al_%i/NMI/k_%i/nmi_max_cons.csv" %>%
              sprintf(path_al, al, k), header = F
          ) %>%
            as.numeric()
          NMI.DF <- NMI.DF %>%
            dplyr::bind_rows(
              dplyr::tibble(
                al = al / 100,
                k = k,
                x = nmi.cos.max,
                type = "NMI"
              )
            )
        }
        if (
          file.exists(
            "%s/al_%i/K/CON/k_%i.csv" %>%
              sprintf(path_al, al, k)
          )
        ) {
          partition.con <- read.csv(
            "%s/al_%i/K/CON/k_%i.csv" %>%
              sprintf(path_al, al, k), header = F
          ) %>%
            t() %>%
            as.numeric()
          PARTITIONS.DF <- PARTITIONS.DF %>%
            dplyr::bind_rows(
              dplyr::tibble(
                alpha = al / 100,
                k = k,
                type = "CON",
                partition = partition.con
              )
            )
        }
        if (
          file.exists(
            "%s/al_%i/K/MAX/k_%i.csv" %>%
              sprintf(path_al, al, k)
          )
        ) {
          partition.max <- read.csv(
            "%s/al_%i/K/MAX/k_%i.csv" %>%
              sprintf(path_al, al, k), header = F
          ) %>%
            t() %>% 
            as.numeric()
          PARTITIONS.DF <- PARTITIONS.DF %>%
            dplyr::bind_rows(
              dplyr::tibble(
                alpha = al / 100,
                k = k,
                type = "MAX",
                partition = partition.max
              )
            )
        }
      }
    }
    GENALX <- rbind(LOGEVCON.DF, LOGEVMAX.DF)
    GENALX$al <- GENALX$al %>%
      as.factor()
    MAX_GENALX <- GENALX %>%
      dplyr::group_by(al, type) %>%
      dplyr::summarise(k = which(x == max(x)) + 1, x = max(x))
    NMI.DF$al <- NMI.DF$al %>%
      as.factor()
    MAX_NMI <- NMI.DF %>%
      dplyr::group_by(al) %>%
      dplyr::summarise(k = which(x == max(x)) + 1, x = max(x))
    # Plots ----
    p <- GENALX %>%
      ggplot2::ggplot(
        ggplot2::aes(color = al, shape = type)
      ) +
      ggplot2::geom_point(
        ggplot2::aes(as.factor(k), x),
        alpha = 0.6,
        size = 3
      ) +
      ggplot2::geom_point(
        data = MAX_GENALX,
        ggplot2::aes(as.factor(k), x, color = al),
        size = 4,
        shape = 6
      ) +
      ggplot2::geom_line(
        ggplot2::aes(k - 1, x),
        alpha = 0.6,
        linetype = "dashed"
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::xlab("k") +
      ggplot2::ylab("LogEvidence") +
      ggplot2::theme_classic()
    q <- ggplot2::ggplot(NMI.DF, ggplot2::aes(color = al)) +
      ggplot2::geom_point(
        ggplot2::aes(as.factor(k), x),
        alpha = 0.6,
        size = 3
      ) +
      ggplot2::geom_point(
        data = MAX_NMI,
        ggplot2::aes(as.factor(k), x, color = al),
        size = 4,
        alpha = 0.7,
        shape = 6
      ) +
      ggplot2::geom_line(
        ggplot2::aes(k - 1, x),
        alpha = 0.6
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::ylab("NMI") +
      ggplot2::xlab("k") +
      ggplot2::theme_classic()
    pp <- ggpubr::ggarrange(
        p, q, labels = "AUTO",
        nrow = 1, ncol = 2
      )
    png(
      "%s/logev_nmi.png" %>%
        sprintf( path.plot),
      width = 12, height = 5, units = "in", res = 200
    )
    print(pp)
    dev.off()
    # Just NMI ----
    q <- ggplot2::ggplot(NMI.DF, ggplot2::aes()) +
      ggplot2::geom_point(
        ggplot2::aes(as.factor(k), x),
        alpha = 0.6,
        size = 3
      ) +
      ggplot2::geom_point(
        data = MAX_NMI,
        ggplot2::aes(as.factor(k), x),
        size = 4,
        alpha = 0.7,
        shape = 6
      ) +
      ggplot2::geom_line(
        ggplot2::aes(k - 1, x),
        alpha = 0.6
      ) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::ylab("NMI") +
      ggplot2::xlab("R") +
      ggplot2::theme_classic()
    png(
      "%s/nmi.png" %>%
        sprintf( path.plot),
      width = 6, height = 5, units = "in", res = 200
    )
    print(q)
    dev.off()
    # Log-EV vs alpha ----
    p <- GENALX %>%
      ggplot2::ggplot(
        ggplot2::aes(color = as.factor(k))
      ) +
      ggplot2::facet_wrap(~type) +
      ggplot2::geom_point(
        ggplot2::aes(al, x),
        alpha = 0.4,
        size = 3
      ) +
      ggplot2::geom_point(
        data = MAX_GENALX,
        ggplot2::aes(al, x, color = as.factor(k)),
        size = 4,
        alpha = 0.8,
        shape = 6
      ) +
      ggplot2::geom_line(
        ggplot2::aes(as.numeric(al), x),
        alpha = 0.6,
        linetype = "dashed"
      ) +
      ggplot2::xlab("alpha") +
      ggplot2::ylab("LogEvidence") +
      ggplot2::ggtitle("LogEvidence vs. al") +
      ggplot2::theme_classic() +
      ggplot2::guides(
        color = ggplot2::guide_legend(title = "k")
      )
    png(
      "%s/logev_alpha.png" %>%
        sprintf(path.plot),
      width = 6, height = 5, units = "in", res = 200
    )
    print(p)
    dev.off()
    NMI.BIG.DF <- data.frame()
    type.1 <- "MAX"
    type.2 <- "CON"
    for (e in K) {
      for (al.1 in alpha / 100) {
        for (al.2 in alpha / 100) {
          if (al.1 <= al.2) {
            part.1 <- PARTITIONS.DF %>%
              dplyr::filter(alpha == al.1, type == type.1, k == e) %>%
              dplyr::select(partition) %>%
              unlist() %>%
              as.numeric()
            part.2 <- PARTITIONS.DF %>%
              dplyr::filter(alpha == al.2, type == type.2, k == e) %>%
              dplyr::select(partition) %>%
              unlist() %>%
              as.numeric()
            var.NMI <- aricode::NMI(part.1, part.2, variant = "sum")
            NMI.BIG.DF <- NMI.BIG.DF %>%
              dplyr::bind_rows(
                dplyr::tibble(
                  alpha.source = al.1,
                  alpha.target = al.2,
                  k = e,
                  type = sprintf("%s-%s", type.1, type.2),
                  nmi = var.NMI
                )
              )
          }
        }
      }
    }
    NMI.BIG.DF$type <- factor(
      NMI.BIG.DF$type,
      levels = c("MAX-MAX", "CON-CON", "MAX-CON")
    )
    NMI.BIG.DF <- NMI.BIG.DF %>%
    dplyr::bind_rows(
      dplyr::tibble(
        alpha.source = NMI.BIG.DF$alpha.target,
        alpha.target = NMI.BIG.DF$alpha.source,
        k = NMI.BIG.DF$k,
        nmi = NMI.BIG.DF$nmi,
        type = NMI.BIG.DF$type
      )
    )
    p <- NMI.BIG.DF %>%
      ggplot2::ggplot(
        ggplot2::aes(
          as.factor(alpha.target),
          as.factor(alpha.source),
          fill = nmi
        )
      ) +
      ggplot2::facet_grid(type ~ k) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis(option = "B") +
      ggplot2::geom_text(
        data = NMI.BIG.DF,
        ggplot2::aes(as.factor(alpha.target), as.factor(alpha.source)),
        label = round(NMI.BIG.DF$nmi,3),
        color = ifelse(NMI.BIG.DF$nmi == 1, "black", "white")
      ) +
      ggplot2::ylab("alpha") +
      ggplot2::xlab("alpha")
    png(
      "%s/alpha_alpha_fill_nmi_facet_type-k.png" %>%
        sprintf(path.plot),
      width = 18, height = 4, units = "in", res = 200
    )
    print(p)
    dev.off()
  } else
    print("No WSBM stats")
}
save.tables.wsbm <- function(K, R, labels, inst, al=0, on=T) {
  if (on) {
    print("Saving WSBM table")
    print("----> Warning: Choose the right partition")
    for (r in R) {
      paste(
        "../CSV", inst$mermoved, "tables", inst$common, K,
        sep = "/"
      ) %>%
        dir.create(showWarnings = F, recursive = T)
      commships <- "../WSBM/CD/CSV/labels/%s/%s/%s/al_0/K/CON/k_%i.csv" %>%
        sprintf(inst$data, inst$mermoved, inst$common, r) %>%
        read.csv(header = F) %>%
        unlist() %>%
        as.numeric()
      dplyr::tibble(
        AREA = labels,
        COMMSHIP = commships,
        DIR = "WSBM"
      ) %>%
        write.csv(
          paste(
            "../CSV",
            inst$mermoved,
            "tables",
            inst$common,
            K,
            "wsbm_al_%i_CON_k_%i.csv" %>%
              sprintf(al, r),
            sep = "/"
          )
        )
    }
  } else{
    print("No WSBM table")
  }
}

main <- function(inst) {
  library(magrittr)
  source("functions/load_net.R")
  netx <- load.net(inst)
  net <- netx$net
  nodes <- netx$nodes
  labels <- netx$labels
  net <- net %>%
    dplyr::filter(source <= nodes)
  source("functions/sformat.R")
  inst$mainfolder <- paste(
    inst$mermoved,
    inst$common,
    paste("WSBM", data, sep = "_"), sep = "/"
  )
  dir.create(
    "%s/%s" %>%
      sprintf(inst$plot, inst$mainfolder),
    showWarnings = F
  )
  path_al <- "../WSBM/CD/CSV/labels/%s/%s/%s" %>%
    sprintf(inst$data, inst$mermoved, inst$common)

  ### Start algorithms
  Data <- make.WSBM.matrix(
    net, labels, nodes, path_al, inst,
    on = F, plt = F
  )
  make.WSBM.stats(path_al, inst, on = T)
  ### Save tables
  save.tables.wsbm(4, 6, labels, inst, al = 0, on = T)
}

data <- "commships"
mermoved <- "merged"
distances <- "tracto2016"
model <- "normal_GB_GB"
csv.name <- "fln_norm"
labels.name <- "imputation_labels"
common.path <- paste(distances, model, sep = "/")
csv.path <- paste(
  mermoved,
  "imputation",
  common.path,
  paste0(csv.name, ".csv"), sep = "/"
)
labels.path <- paste(
  mermoved,
  "labels",
  common.path,
  paste0(labels.name, ".csv"), sep = "/"
)
plot.path <- "../plots"

path.list <- list(
  csv = csv.path,
  labels = labels.path,
  plot = plot.path,
  common = common.path,
  data = data,
  mermoved = mermoved,
  model = model,
  distances = distances
)

main(path.list)

### Code to plot networks with points and compare the max and consensus plots
### using the Munkres algorithm
# source("funtions/apply_munkres.R")
# for (sub in c('MAX', 'CON')){
#   for (k in 2:7){
#     reference.membership <- Data$membership[Data$k == k &
#                                               Data$alpha == 0.75 &
#                                               Data$type == 'MAX']
#     for (al in c(0, 0.25, 0.5, 0.75, 1)){
#       old.memberships <- Data$membership[Data$k == k &
#                                          Data$alpha == al &
#                                          Data$type == sub]
# 
#       Data$membership[Data$k == k & Data$alpha == al & Data$type == sub] <- apply.munkres(old.memberships, reference.membership, k)
#     }
#   }
# }
# 
# Data$membership <- Data$membership %>% as.factor()
# Data$k <- Data$k %>% as.factor()
# Data$alpha <- Data$alpha %>% as.factor()
# 
# p <- ggplot(Data %>% filter(type=='MAX'), aes(x,y, fill=membership))+
#   facet_grid(k ~ alpha)+
#   geom_point(shape=21, size=7, stroke = 0.8)+
#   geom_text(label=Data$labels[Data$type == 'MAX'], color='white', size=1.5, fontface='bold') +
#   ggtitle('MAX logEV')+
# ggplot(Data %>% filter(type=='CON'), aes(x,y, fill=membership))+
#   facet_grid(k ~ alpha)+
#   geom_point(shape=21, size=7, stroke = 0.8)+
#   geom_text(label=Data$labels[Data$type == 'CON'], color='white', size=1.5, fontface='bold')+
#   ggtitle('CON logEv')
# 
# cairo_ps(filename = sprintf("/Volumes/JMA_1/FLN_52_108/plots/WSBM/Summaries/%s.eps", 'CON_MAX_clustering_all_k'),
#          width =40, height = 20,
#          fallback_resolution = 200)
# print(p)
# dev.off()
# 
# p <- ggplot(Data %>% filter(type=='MAX'), aes(x,y, fill=membership))+
#   facet_grid(k ~ alpha)+
#   geom_point(shape=21, size=7, stroke = 0.8)+
#   geom_text(label=Data$labels[Data$type == 'MAX'], color='white', size=1.5, fontface='bold') +
#   ggtitle('MAX logEV')
# 
# cairo_ps(filename = sprintf("/Volumes/JMA_1/FLN_52_108/plots/WSBM/Summaries/%s.eps", 'MAX_clustering_all_k'),
#          width =24, height = 20,
#          fallback_resolution = 200)
# print(p)
# dev.off()
# 
# p <- ggplot(Data %>% filter(type=='CON'), aes(x,y, fill=membership))+
#   facet_grid(k ~ alpha)+
#   geom_point(shape=21, size=7, stroke = 0.8)+
#   geom_text(label=Data$labels[Data$type == 'CON'], color='white', size=1.5, fontface='bold') +
#   ggtitle('MAX logEV')
# 
# cairo_ps(filename = sprintf("/Volumes/JMA_1/FLN_52_108/plots/WSBM/Summaries/%s.eps", 'CON_clustering_all_k'),
#          width =24, height = 20,
#          fallback_resolution = 200)
# print(p)
# dev.off()
