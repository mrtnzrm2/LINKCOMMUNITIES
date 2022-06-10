standardized <- function(v) {
  return((v - mean(v, na.rm = T)) / sd(v, na.rm = T))
}
standardized.fln <- function(v) {
  v[v == 0] <- NA
  return((v - mean(v, na.rm = T)) / sd(v, na.rm = T))
}

standardized.diag <- function(v, nt) {
  v[diag(nt) == 1] <- NA
  return((v - mean(v, na.rm = T)) / sd(v, na.rm = T))
}

standardized.block.diag <- function(v) {
  nt <- dim(v)[2]
  v[1:nt, 1:nt][diag(nt) == 1] <- NA
  return((v - mean(v, na.rm = T)) / sd(v, na.rm = T))
}

divide_sd <- function(v, nt) {
  v[diag(nt) == 1] <- NA
  return(v / sd(v, na.rm = T))
}

divide_max <- function(v, nt) {
  v[diag(nt) == 1] <- NA
  v <- v / max(v, na.rm = T)
  v[diag(nt) == 1] <- 0
  return(v)
}