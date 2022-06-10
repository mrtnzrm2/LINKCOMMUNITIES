plot.heuristics <- function(
  path.plot, foldername, hclust_features, subfolder=""
  ) {
  source("functions/find_height.R")
  hclust_features <- hclust_features[!is.na(hclust_features$SS), ]
  dir.create(
    "%s/%s/%s" %>%
      sprintf(path.plot, foldername, subfolder),
    showWarnings = F
  )
  wrap_hclust <- dplyr::tibble(
    K = hclust_features$K,
    var = hclust_features$Dc,
    val = "Dc"
  ) %>%
    dplyr::bind_rows(
      dplyr::tibble(
        K = hclust_features$K,
        var = hclust_features$SS,
        val = "SS"
      )
    )
  fig <- plotly::plot_ly(data = wrap_hclust, x = ~K) %>%
    plotly::add_trace(
      y = ~var,
      name = ~val,
      mode = "lines+markers",
      color = ~val,
      marker = list(
        size = 3,
        line = list(
          width = 2
        )
      )
    )
  fig <- plotly::layout(fig,  xaxis = list(type = "log"))
  plotly::as_widget(fig) %>%
    htmlwidgets::saveWidget(
      "%s/%s/%s/%s.html" %>%
        sprintf(path.plot, foldername, subfolder, "var_logK")
    )
}