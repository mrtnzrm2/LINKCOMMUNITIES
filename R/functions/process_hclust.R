calculate.s.s <- function(s, M, n=2) {
  s.s <- 0
  ns <- length(s)
  if (ns > 1) {
    if (ns > n) {
      for (i in 2:n) {
        s.s <- s.s + (s[i] * (s[i] + s[i - 1])) / (s[i - 1] * M)
      }
      s.s <- s.s / (n - 1)
    } else{
      for (i in 2:ns) {
        s.s <- s.s + (s[i] * (s[i] + s[i - 1])) / (s[i - 1] * M)
      }
      s.s <- s.s / (ns - 1)
    }
  } else
    s.s <- NA
  return(s.s)
}

solve_M <- function(M) {
  x1 <- (1 + sqrt(1 + (4 * M))) / 2 %>%
    floor()
  x2 <- (1 - sqrt(1 + (4 * M))) / 2 %>%
    floor()
  if (x2 <= 0)
    x <- x1
  else
    x <- x2
  return(x)
}

process.hclust <- function(net, net.cluster,  nodes, M, type = "directed") {
  height <- net.cluster$height[!duplicated(net.cluster$height)]
  leaves <- length(height)
  Dc.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  Nmin.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  Nmax.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  K.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  NEC.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  NAC.view <- matrix(nodes, nrow = leaves + 1, ncol = 1)
  s.s.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  K.view[1] <- max(net$id)

  for (i in 1:leaves) {
    cluster.h <- cutree(net.cluster, h = height[i])
    clusters <- unique(cluster.h)
    K.view[i + 1] <- length(clusters)
    Dc <- matrix(0, nrow = length(clusters), ncol = 1)
    Mc <- matrix(0, nrow = length(clusters), ncol = 1)
    Nc <- matrix(2, nrow = length(clusters), ncol = 1)
    NEC <- matrix(0, nrow = length(clusters), ncol = 1)
    for (ii in 1:length(clusters)) {
      net.cls <- net[cluster.h == clusters[ii], c("source", "target")]
      Mc[ii] <- nrow(net.cls)
      sor <- unique(net.cls$source)
      tar <- unique(net.cls$target)
      ndc <- length(unique(c(sor, tar)))
      if (ndc > 2 && Mc[ii] > 1) {
        if (type == "directed") {
          Dc[ii] <- (Mc[ii] - ndc + 1) / ((ndc - 1)**2)
          NEC[ii] <- 1
          Nc[ii] <- ndc
          NAC.view[i + 1] <- NAC.view[i + 1] - solve_M(Mc[ii])
        } else if (type == "undirected") {
          Dc[ii] <- 2 * (Mc[ii] - ndc + 1) / ((ndc - 1) * (ndc - 2))
          NEC[ii] <- 1
          Nc[ii] <- ndc
          NAC.view[i + 1] <- NAC.view[i + 1] - solve_M(2 * Mc[ii])
        }
      }
    }
    s.s <- cluster.h %>%
      table() %>%
      sort(decreasing = T)
    s.s <- s.s %>%
      unname()
    if (i > leaves - 10)
      print(s.s)
    Dc.view[i + 1] <- t(Mc) %*% Dc / M
    Nmin.view[i + 1] <- min(Nc, na.rm = T)
    Nmax.view[i + 1] <- max(Nc, na.rm = T)
    NEC.view[i + 1] <- sum(NEC == 1)
    s.s.view[i + 1] <- calculate.s.s(s.s, M, n = 6)
  }
  height <- c(0, height)
  D.network <- dplyr::tibble(
    "height" = height,
    "Dc" = Dc.view,
    "Nmin" = Nmin.view,
    "Nmax" = Nmax.view,
    "K" = K.view,
    "NEC" = NEC.view,
    "NAC" = NAC.view,
    "SS" = s.s.view
  )
  return(D.network)
}