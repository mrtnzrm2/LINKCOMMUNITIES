format.regions <- function(regions){
  for (i in 1:nrow(regions)){
    if (grepl("-", regions$AREA[i], fixed = T)){
      regions$AREA[i] <- stringr::str_replace(regions$AREA[i], "-", ".")
    }
    if (grepl(" ", regions$AREA[i], fixed = T)){
      regions$AREA[i] <- stringr::str_replace(regions$AREA[i], " ", "")
    }
    if (regions$AREA[i] == "aud.core")
      regions$AREA[i] <- "core"
    if (regions$AREA[i] == "opal")
      regions$AREA[i] <- "opai"
    if (regions$AREA[i] == "t.pole")
      regions$AREA[i] <- "pole"
    if (regions$AREA[i] == "prost")
      regions$AREA[i] <- "pro.st."
    if (regions$AREA[i] == "35.36")
      regions$AREA[i] <- "peri"
  }
  regions <- regions %>% rbind(data.frame(AREA="pi", REGION="Temporal"))
  return(regions)
}