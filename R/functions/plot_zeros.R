plot.zeros <- function(train, test, kcoords, serie, inst, subfolder="") {
    dir.create(
        "%s/%s/Regression/XGBOOST/%s/%i" %>%
            sprintf(inst$plot, inst$folder, subfolder, serie),
        showWarnings = F
    )
     # Joining zeros informations ----
    train_link <- train %>%
    dplyr::rename(W = w)
    train_link$cat <- ifelse(train_link$W > 0, "1", "0")
    train_link$cat.2 <- "Train"
    test_link <- test %>%
    dplyr::rename(W = w)
    test_link$cat <- ifelse(test_link$W > 0, "1", "0")
    test_link$cat.2 <- "Test"
    data <- train_link %>%
    dplyr::bind_rows(test_link)
    kcoords <- kcoords %>%
        dplyr::as_tibble() %>%
        dplyr::rename(
            dist = V2,
            sim = V1
        )
    # Scatter-plot dist-sim
    p_scat <- data %>%
        ggplot2::ggplot(ggplot2::aes(dist, sim, color = cat)) +
        ggplot2::facet_wrap(~cat.2) +
        ggplot2::geom_point(size = 0.5, alpha = 0.5) +
        ggplot2::geom_point(
            data = kcoords,
            ggplot2::aes(dist, sim),
            color = "black"
        ) +
        ggplot2::theme_bw()
    # Histogram sim
    p_histosim <- data %>%
        ggplot2::ggplot(ggplot2::aes(sim, fill = cat)) +
        ggplot2::facet_wrap(~cat.2) +
        ggplot2::geom_histogram(
            ggplot2::aes(y = ..density..),
            position = "identity",
            color = "gray",
            alpha = 0.5,
        ) +
        ggplot2::theme_bw()
    # Histogram dist
    p_histodist <- data %>%
        ggplot2::ggplot(ggplot2::aes(dist, fill = cat)) +
        ggplot2::facet_wrap(~cat.2) +
        ggplot2::geom_histogram(
            ggplot2::aes(y = ..density..),
            position = "identity",
            color = "gray",
            alpha = 0.5,
        ) +
        ggplot2::theme_bw()
    p <- cowplot::plot_grid(p_scat, p_histosim, p_histodist, nrow = 3)
    png("%s/%s/Regression/XGBOOST/%s/%i/zero_information.png" %>%
        sprintf(inst$plot, inst$folder, subfolder, serie),
        width = 10, height = 15, res = 200, units = "in"
    )
    print(p)
    dev.off()
}