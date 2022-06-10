arrange.memberships.values <- function(memberships) {
  umembs <- memberships %>%
    unique() %>%
    sort()
  rmembs <- seq_len(length(umembs))
  kmax <- max(memberships)
  if (length(umembs) != kmax)
    print("** Warning: There are less node-communities tham the maximum K")
  for (i in seq_len(length(memberships))) {
    memberships[i] <- rmembs[umembs == memberships[i]]
  }
  return(memberships)
}