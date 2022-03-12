get.bourder.regions <- function(regions){
  
  nr <- regions$REGION %>% unique() %>% length()
  bourder <- rep("", nr-1)
  current <- regions$REGION[1]
  count <- 1
  for (i in 2:nrow(regions)){
    if(current != regions$REGION[i]){
      bourder[count] <- regions$AREA[i]
      current <- regions$REGION[i]
      count <- count + 1
    }
  }
  return(bourder)
}