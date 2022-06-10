arrange.memberships.reference <- function(
  memberships, k, inst
) {
  ref.membership <- "../WSBM/CD/CSV/labels/%s/%s/%s/al_%i/K/%s/k_%i.csv" %>%
    sprintf(
      inst$data, inst$folder, inst$common,
      inst$al, inst$MAN, k
    ) %>%
    read.csv() %>%
    unlist() %>%
    as.numeric()
  source("functions/apply_munkres.R")
  memberships <- apply.munkres(memberships, ref.membership, k)
  return(memberships)
}