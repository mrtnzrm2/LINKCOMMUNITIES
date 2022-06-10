get.bourder.labels <- function(memberships, labels) {
  get.order <- order(memberships)
  memberships <- memberships[get.order]
  labels <- labels[get.order]
  nm <- memberships %>%
    unique() %>%
    length()
  tab.mem <- memberships %>%
    table()
  bourders <- rep("", nm - 1)
  for (i in 1:nm) {
    bourders[i] <- labels[sum(tab.mem[1:i])]
  }
  return(bourders[1:(length(bourders) - 1)])
}