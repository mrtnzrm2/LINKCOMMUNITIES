library(magrittr)
model <- "normal_GB_GB"
fln <- read.csv(
  "../CSV/merged/imputation/tracto2016/%s/fln_raw.csv" %>%
    sprintf(model)
)
adj <- read.csv(
  "../CSV/merged/imputation/tracto2016/%s/adj_raw.csv" %>%
    sprintf(model)
)
fln_colnames <- fln %>%
  colnames()
adj_colnames <- adj %>%
  colnames()
adj_colnames <- adj_colnames[1:107]
adj <- adj[, 1:107]
neworder_adj <- match(fln_colnames, adj_colnames)
fln <- fln %>%
  as.matrix()
adj <- adj %>%
  as.matrix()
adj <- adj[, neworder_adj]
adj <- adj[neworder_adj, ]
adj[adj <= 0.418138] <- 0

fln_adj <- fln * adj
fln_adj_norm <- sweep(
  fln_adj,
  2,
  apply(fln_adj, 2, sum),
  FUN = "/"
)

write.csv(
  fln_adj,
  "../CSV/merged/imputation/tracto2016/%s/fln.csv" %>%
    sprintf(model),
  row.names = F
)
write.csv(
  fln_adj_norm,
  "../CSV/merged/imputation/tracto2016/%s/fln_norm.csv" %>%
    sprintf(model),
  row.names = F
)
# source("functions/adj_to_df.R")
# fln_adj <- fln_adj %>%
#   adj.to.df()
# fln_adj <- fln_adj[fln_adj$weight > 0, ]
# fln_adj$type <- ifelse(
#   fln_adj$target <= 50,
#   "INJ",
#   "MISS"
# )
# mean_fln_adj <- fln_adj %>%
#   dplyr::group_by(type) %>%
#   dplyr::summarise(
#     mean = mean(weight) %>%
#       round(3)
#   )
# ####
# fln_adj_norm <- fln_adj_norm %>%
#   adj.to.df()
# fln_adj_norm <- fln_adj_norm[
#   fln_adj_norm$weight > 0,
# ]
# fln_adj_norm$type <- ifelse(
#   fln_adj_norm$target <= 50,
#   "INJ",
#   "MISS"
# )
# mean_fln_adj_norm <- fln_adj_norm %>%
#   dplyr::group_by(type) %>%
#   dplyr::summarise(
#     mean = mean(weight) %>%
#       round(3)
  # )
# Plots ----
# p <- fln_adj %>%
#   ggplot2::ggplot(
#     ggplot2::aes(weight, fill = type)
#   ) +
#   ggplot2::geom_histogram(
#     ggplot2::aes(
#       y = ..density..
#     ),
#     position = "identity",
#     color = "black",
#     alpha = 0.4
#   ) +
#   ggplot2::annotation_custom(
#     gridExtra::tableGrob(mean_fln_adj),
#     ymin = 0.3
#   ) +
#   ggplot2::ggtitle("Unnormalized") +
#   ggplot2::scale_x_continuous(trans = "log10") +
#   ggplot2::theme_classic() +
#   ggplot2::theme(
#     legend.position = "none"
#   )
# p_norm <- fln_adj_norm %>%
#   ggplot2::ggplot(
#     ggplot2::aes(weight, fill = type)
#   ) +
#   ggplot2::geom_histogram(
#     ggplot2::aes(
#       y = ..density..
#     ),
#     position = "identity",
#     color = "black",
#     alpha = 0.4
#   ) +
#   ggplot2::annotation_custom(
#     gridExtra::tableGrob(mean_fln_adj_norm),
#     xmin = 0.02,
#     ymin = 0.3
#   ) +
#   ggplot2::ggtitle("Normalized") +
#   ggplot2::scale_x_continuous(trans = "log10") +
#   ggplot2::theme_classic()
# P <- cowplot::plot_grid(
#   p, p_norm, ncol = 2, rel_widths = c(1, 1.3), labels = "AUTO"
# )
# png(
#   "../plots/merged/PreAnalysis/normalized_or_unnormalized.png",
#   width = 14, height = 5, units = "in", res = 200
# )
# print(P)
# dev.off()
# fln_org <- read.csv(
#   "../CSV/merged/CountMatrix_Summed_Main43_Periphery8_220126_arr_fln.csv"
# ) %>%
#   as.matrix()
# fln_org_colnames <- fln_org %>%
#   colnames() %>%
#   tolower()
# fln_colnames <- fln_colnames %>%
#   tolower()
# print(
#   fln_org_colnames[!fln_org_colnames %in% fln_colnames[1:50]]
# )
# print(
#   fln_colnames[1:50][!fln_colnames[1:50] %in% fln_org_colnames]
# )