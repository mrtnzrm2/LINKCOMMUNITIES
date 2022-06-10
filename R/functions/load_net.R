format_areas_regions <- function(areas_region) {
  areas_region[areas_region == "29-30"] <- "29/30"
  areas_region[areas_region == "prost"] <- "pro.st."
  areas_region[areas_region == "tea-m a"] <- "tea/ma"
  areas_region[areas_region == "tea-m p"] <- "tea/mp"
  areas_region[areas_region == "th-tf"] <- "th/tf"
  areas_region[areas_region == "9-46d"] <- "9/46d"
  areas_region[areas_region == "9-46v"] <- "9/46v"
  areas_region[areas_region == "opal"] <- "opai"
  areas_region[areas_region == "t.pole"] <- "pole"
  areas_region[areas_region == "parains"] <- "pi"
  areas_region[areas_region == "insula"] <- "ins"
  areas_region[areas_region == "aud. core"] <- "core"
  areas_region[areas_region == "35-36"] <- "35/36"
  return(areas_region)
}

format_areas_matrix <- function(areas_matrix) {
  areas_matrix[areas_matrix == "9.46d"] <- "9/46d"
  areas_matrix[areas_matrix == "9.46v"] <- "9/46v"
  areas_matrix[areas_matrix == "29.30"] <- "29/30"
  areas_matrix[areas_matrix == "insula"] <- "ins"
  areas_matrix[areas_matrix == "tea.ma"] <- "tea/ma"
  areas_matrix[areas_matrix == "tea.mp"] <- "tea/mp"
  areas_matrix[areas_matrix == "th.tf"] <- "th/tf"
  return(areas_matrix)
}

regions_lobes <- function(labels) {
  source("functions/gg_color_hue.R")
  regions <- read.csv(
    file.path("..", "CSV/Regions/Table_areas_regions_09_2019.csv")
  ) %>%
    dplyr::as_tibble()
  colnames(regions) <- c("AREA", "REGION")
  regions$AREA <- regions$AREA %>%
    tolower() %>%
    format_areas_regions()
  regions <- regions[regions$AREA %in% labels, ]
  lobes <- regions$REGION %>%
    unique() %>%
    sort()
  color_region <- dplyr::tibble(
    REGION = c(
      "Occipital",
      "Temporal",
      "Parietal",
      "Frontal",
      "Prefrontal",
      "Cingulate"
    ),
    COLOR = c(
      rgb(0, 97, 62, maxColorValue = 255),
      rgb(255, 126, 0, maxColorValue = 255),
      "#800080",
      "#FFD500",
      rgb(237, 28, 36, maxColorValue = 255),
      "#2A52BE"
    )
  )
  color_region <- color_region[order(color_region$REGION), ]
  regions$COLOR <- color_region$COLOR[
    match(
      regions$REGION, color_region$REGION
    )
  ]
  return(regions)
}

regions_streams <- function() {
  dorsal <- c(
    "7A", "STPc", "STPi", "STPr", "FST", "MST",
    "MT", "LIP", "PIP", "MIP", "VIP", "PGa",
    "V3A", "DP", "V6", "V6A", "V4pcLF", "V4fpLF"
  ) %>%
    tolower()
  ventral <- c(
    "TH/TF", "TEav", "TEad", "TEa/mp",
    "TEa/ma", "TEpv", "TEpd", "TEOm", "TEO",
    "PERI", "V4c", "V4pcUF", "V4fpUF"
  ) %>%
    tolower()
  sudo.dorsal <- c("V1pcLF", "V1fpLF") %>%
    tolower()
  sudo.ventral <- c("V1fpUF", "V2fpUF", "V2pcUF", "V1c") %>%
    tolower()
  regions <- data.frame(AREA = labels, REGION = NA)
  regions$REGION[labels %in% dorsal] <- "dorsal"
  regions$REGION[labels %in% ventral] <- "ventral"
  regions$REGION[labels %in% sudo.dorsal] <- "sudo_dorsal"
  regions$REGION[labels %in% sudo.ventral] <- "sudo_ventral"
  regions$REGION[is.na(regions$REGION)] <- "other"
  regions$COLOR <- NA
  regions$COLOR[regions$REGION == "dorsal"] <-
    rgb(1, 0.2, 0, 0.9)
  regions$COLOR[regions$REGION == "ventral"] <-
    rgb(25, 25, 112, 230, maxColorValue = 255)
  regions$COLOR[regions$REGION == "sudo_dorsal"] <-
    rgb(233, 150, 122, 230, maxColorValue = 255)
  regions$COLOR[regions$REGION == "sudo_ventral"] <-
    rgb(100, 149, 237, 230, maxColorValue = 255)
  regions$COLOR[regions$REGION == "other"] <-
    rgb(60, 179, 113, 230, maxColorValue = 255)
  return(regions)
}

load.net <- function(inst) {
  library(magrittr)
  source("functions/sformat.R")
  net.fln <- read.csv(file.path("..", "CSV", inst$csv))
  fln.coords<- read.csv(file.path("..", "CSV", "fln_2d.csv")) %>%
    as.matrix()
  dis.coords<- read.csv(file.path("..", "CSV", "dis_2d.csv")) %>%
    as.matrix()
  net.fln <- as.matrix(net.fln)
  net.fln[net.fln == 0] <- NA
  source("functions/format_labels.R")
  labels <- read.csv(file.path("..", "CSV", inst$labels)) %>%
    as.matrix() %>%
    unname() %>%
    tolower() %>%
    format.labels()
  if (grepl("merged", inst$csv, fixed = T))
    labels <- labels %>%
      format_areas_matrix()
  leaves <- length(net.fln[!is.na(net.fln)])
  source("functions/adj_to_df.R")
  net <- net.fln %>%
    adj.to.df()
  net <- net[!is.na(net$weight), ]
  net$id <- 1:leaves
  nodes <- max(net$target)
  net$fln <- net$weight
  net$w <- log10(net$weight) + 7
  # regions <- regions_streams()
  regions <- regions_lobes(labels)

  x <- list()
  x$net <- net
  x$nodes <- nodes
  # x$leaves <- length(net$source[net$source <= nodes])
  x$leaves <- leaves
  x$labels <- labels
  x$fln.coords <- fln.coords
  x$dis.coords <- dis.coords
  x$regions <- regions
  return(x)
}