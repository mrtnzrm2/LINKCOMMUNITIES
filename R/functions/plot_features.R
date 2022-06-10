plot.features <- function(nt, names, ..., path="", folder="", subfolder="") {
  # All inputs as matrices
  # Organize data ----
  arguments <- list(...)
  na <- length(arguments)
  for (i in seq_len(na)) {
    arguments[[i]] <- arguments[[i]][1:nt, 1:nt] %>%
    t() %>%
    pracma::Reshape(nt * nt, 1)
  }
  data <- arguments[[1]] %>%
    dplyr::as_tibble()
  for (i in 2:na) {
    data <- data %>%
      dplyr::bind_cols(
        arguments[[i]] %>%
          dplyr::as_tibble()
      )
  }
  colnames(data) <- names
  print(min(data$jaccp_src, na.rm = T))
  print(min(data$jaccp_tgt, na.rm = T))
  # Scatter plots ----
  ## Scatter jaccp_tgt vs. jacc_src ----
  scat.src.tgt <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(jaccp_src, jaccp_tgt, color = w_standard)
    ) +
      ggplot2::geom_point(size = 0.5) +
      viridis::scale_color_viridis(na.value = "grey") +
      ggplot2::geom_smooth(method = "lm", color = "red") +
      ggplot2::xlab("standardized jaccard source") +
      ggplot2::ylab("standardized jaccard target") +
      ggpubr::stat_cor() +
      ggplot2::theme_bw()
  ## Scatter dist vs. jaccp_src ----
  scat.src.dist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(jaccp_src, dist, color = w_standard)
    ) +
      ggplot2::geom_point(size = 0.5) +
      viridis::scale_color_viridis(na.value = "grey") +
      ggplot2::geom_smooth(method = "lm", color = "red") +
      ggplot2::xlab("standardized jaccard source") +
      ggplot2::ylab("standardized distance") +
      ggpubr::stat_cor() +
      ggplot2::theme_bw()
  ## Scatter dist vs. jaccp_src ----
  scat.tgt.dist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(jaccp_tgt, dist, color = w_standard)
    ) +
      ggplot2::geom_point(size = 0.5) +
      viridis::scale_color_viridis(na.value = "grey") +
      ggplot2::geom_smooth(method = "lm", color = "red") +
      ggplot2::xlab("standardized jaccard target") +
      ggplot2::ylab("standardized distance") +
      ggpubr::stat_cor() +
      ggplot2::theme_bw()
  ## Scatter dist vs. w, color = jaccp_src ----
  scat.dist.w.src <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, w_standard, color = jaccp_src)
    ) +
      ggplot2::geom_point(size = 0.5) +
      viridis::scale_color_viridis(na.value = "grey") +
      ggplot2::geom_smooth(method = "lm", color = "red") +
      ggplot2::xlab("standardized distance") +
      ggplot2::ylab("Weight") +
      ggpubr::stat_cor() +
      ggplot2::theme_bw()
  ## Scatter dist vs. w, color = jaccp_tgt ----
  scat.dist.w.tgt <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, w_standard, color = jaccp_tgt)
    ) +
      ggplot2::geom_point(size = 0.5) +
      viridis::scale_color_viridis(na.value = "grey") +
      ggplot2::geom_smooth(method = "lm", color = "red") +
      ggplot2::xlab("standardized distance") +
      ggplot2::ylab("Weight") +
      ggpubr::stat_cor() +
      ggplot2::theme_bw()
  ## Scatter jaccp_src vs. w, color = distance ----
  scat.src.w.dist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(jaccp_src, w_standard, color = jaccp_tgt)
    ) +
      ggplot2::geom_point(size = 0.5) +
      viridis::scale_color_viridis(na.value = "grey") +
      ggplot2::geom_smooth(method = "lm", color = "red") +
      ggplot2::xlab("standardized jaccard source") +
      ggplot2::ylab("Weight") +
      ggpubr::stat_cor() +
      ggplot2::theme_bw()
    p.scatter <- cowplot::plot_grid(
    scat.src.tgt, scat.src.dist, scat.tgt.dist,
    scat.dist.w.src, scat.dist.w.tgt, scat.src.w.dist,
    ncol = 3, nrow = 2, labels = "auto"
  )
  # Histograms ----
  hist.similarity <- dplyr::tibble(
    average_jaccp = rowMeans(
      cbind(data$jaccp_src, data$jaccp_tgt), na.rm = T
    ),
    exist = ifelse(is.na(data$w_standard), "0", "1")
  )
  mean.hist.sim <- hist.similarity %>%
    dplyr::group_by(exist) %>%
    dplyr::summarise(
      average_jaccp = mean(average_jaccp, na.rm = T)
    )
  hist.src <- dplyr::tibble(
    jaccp_src = data$jaccp_src,
    exist = ifelse(is.na(data$w_standard), "0", "1")
  )
  mean.hist.src <- hist.src %>%
    dplyr::group_by(exist) %>%
    dplyr::summarise(
      jaccp_src = mean(jaccp_src, na.rm = T)
    )
  hist.tgt <- dplyr::tibble(
    jaccp_tgt = data$jaccp_tgt,
    exist = ifelse(is.na(data$w_standard), "0", "1")
  )
  mean.hist.tgt <- hist.tgt %>%
    dplyr::group_by(exist) %>%
    dplyr::summarise(
      jaccp_tgt = mean(jaccp_tgt, na.rm = T)
    )
  hist.dist <- dplyr::tibble(
    dist = data$dist,
    exist = ifelse(is.na(data$w_standard), "0", "1")
  )
  mean.hist.dist <- hist.dist %>%
    dplyr::group_by(exist) %>%
    dplyr::summarise(
      dist = mean(dist, na.rm = T)
    )
  hist.w <- dplyr::tibble(
    w = data$w_standard
  )
  hist.similarity <- hist.similarity %>%
    ggplot2::ggplot(
      ggplot2::aes(average_jaccp, fill = exist)
    ) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ..density..),
        position = "identity", bins = 50, alpha = 0.5, color = "gray"
      ) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::geom_vline(
        xintercept = mean.hist.sim$average_jaccp,
        alpha = 0.5,
        linetype = "dashed"
      ) +
      ggplot2::geom_text(
        data = mean.hist.sim,
        ggplot2::aes(average_jaccp + 0.5, 0.65),
        label = mean.hist.sim$average_jaccp %>%
          round(2)
      ) +
      ggplot2::xlab("standardized average jaccard") +
      ggplot2::theme_bw()
  hist.src <- hist.src %>%
    ggplot2::ggplot(
      ggplot2::aes(jaccp_src, fill = exist)
    ) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ..density..),
        position = "identity", bins = 50, alpha = 0.5, color = "gray"
      ) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::geom_vline(
        xintercept = mean.hist.src$jaccp_src,
        alpha = 0.5,
        linetype = "dashed"
      ) +
      ggplot2::geom_text(
        data = mean.hist.src,
        ggplot2::aes(jaccp_src + 0.5, 0.65),
        label = mean.hist.src$jaccp_src %>%
          round(2)
      ) +
      ggplot2::xlab("standardized jaccard source") +
      ggplot2::theme_bw()
  hist.tgt <- hist.tgt %>%
    ggplot2::ggplot(
      ggplot2::aes(jaccp_tgt, fill = exist)
    ) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ..density..),
        position = "identity", bins = 50, alpha = 0.5, color = "gray"
      ) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::geom_vline(
        data = mean.hist.tgt,
        xintercept = mean.hist.tgt$jaccp_tgt,
        alpha = 0.5,
        linetype = "dashed"
      ) +
      ggplot2::geom_text(
        data = mean.hist.tgt,
        ggplot2::aes(jaccp_tgt + 0.5, 0.65),
        label = mean.hist.tgt$jaccp_tgt %>%
          round(2)
      ) +
      ggplot2::xlab("standardized jaccard target") +
      ggplot2::theme_bw()
  hist.dist <- hist.dist %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, fill = exist)
    ) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ..density..),
        position = "identity", bins = 30, alpha = 0.5, color = "gray"
      ) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::geom_vline(
        xintercept = mean.hist.dist$dist,
        alpha = 0.5,
        linetype = "dashed"
      ) +
      ggplot2::geom_text(
        data = mean.hist.dist,
        ggplot2::aes(dist + 0.5, 0.65),
        label = mean.hist.dist$dist %>%
          round(2)
      ) +
      ggplot2::xlab("standardized distance") +
      ggplot2::theme_bw()
  hist.w <- hist.w %>%
    ggplot2::ggplot(
      ggplot2::aes(w)
    ) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ..density..),
        position = "identity", bins = 30,
        alpha = 0.5, color = "gray", fill = "blue"
      ) +
      ggplot2::geom_density(alpha = 0.3, fill = "blue") +
      ggplot2::stat_function(fun = dnorm, color = "red") +
      ggplot2::xlab("standardized weight") +
      ggplot2::theme_bw()
  p.hist <- cowplot::plot_grid(
    hist.similarity, hist.src, hist.tgt, hist.dist, hist.w,
    nrow = 2, ncol = 3, labels = "auto"
  )
  # Plot ----
  png(
    "%s/%s/%s/scatter_features.png" %>%
      sprintf(path, folder, subfolder),
    width = 15, height = 10, units = "in", res = 200
  )
  print(p.scatter)
  dev.off()
  png(
    "%s/%s/%s/histogram_features.png" %>%
      sprintf(path, folder, subfolder),
    width = 15, height = 10, units = "in", res = 200
  )
  print(p.hist)
  dev.off()
}