get.regions.colors <- function(labels, nt){
  ## Get and format regions
  print("** Warning: Be careful selecting the right regions")
  regions <- read.csv("../CSV/Regions/Table_areas_regions_09_2019.csv")
  colnames(regions) <- c("AREA", "REGION")
  regions$AREA <- regions$AREA %>% tolower()
  
  source("functions/format_regions.R")
  regions <- format.regions(regions)
  regions <- regions[order(regions$REGION),]
  
  source("functions/get_bourder_regions.R")
  bourder.out <- get.bourder.regions(regions)
  bourder.in <- get.bourder.regions(regions[regions$AREA %in% labels[1:nt],])
  
  ## Set color for regions
  source("functions/gg_color_hue.R")
  table.regions.out <- table(regions$REGION[regions$AREA %in% labels])
  nr <- table.regions.out %>% names() %>% length()
  color.regions <- gg.color.hue(nr)
  colors.out <- c()
  for (i in 1:nr){
    colors.out <- c(colors.out, rep(color.regions[i], table.regions.out[i]))
  }
  table.regions.in <- table(regions$REGION[regions$AREA %in% labels[1:nt]])
  nr <- table.regions.in %>% names() %>% length()
  color.regions <- gg.color.hue(nr)
  colors.in <- c()
  for (i in 1:nr){
    colors.in <- c(colors.in, rep(color.regions[i], table.regions.in[i]))
  }
  
  x <- list(cin=colors.in, cout=colors.out, bin=bourder.in, bout=bourder.out, rgn=regions)
  return(x)
}