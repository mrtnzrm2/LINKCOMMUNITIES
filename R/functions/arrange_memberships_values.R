arrange.memberships.values <- function(memberships){
  
  umembs <- memberships %>% unique() %>% sort()
  rmembs <- 1:length(umembs)
  kmax <- max(memberships)
  
  if (length(umembs) != kmax)
    print("** Warning: There are less node-communities tham the maximum K")
  
  for (i in 1:length(memberships)){
    memberships[i] <- rmembs[umembs == memberships[i]]
  }
  
  
  return(memberships)
}