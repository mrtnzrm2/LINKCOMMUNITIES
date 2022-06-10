split.dataset.kfold <- function(N, kfold=5) {
  nodes <- 1:N
  ntrn <- N * ((kfold - 1) / kfold)
  ntrn <- ntrn %>%
    as.integer()
  ntst <- N - ntrn
  node.resample <- sample(nodes, N, replace = F)
  train.tb <- dplyr::tibble()
  test.tb <- dplyr::tibble()
  for (i in 1:kfold) {
    x <- (1 + ntst * (i - 1))
    y <- (ntst * i)
    if (i == kfold && y > N)
      y <- N
    test.nodes <- node.resample[x:y]
    train.nodes <- nodes[!nodes %in% test.nodes]
    train.tb <- train.tb %>%
      dplyr::bind_rows(
        dplyr::tibble(nodes = train.nodes, fold = i)
      )
    test.tb <- test.tb %>%
      dplyr::bind_rows(
        dplyr::tibble(nodes = test.nodes, fold = i)
      )
  }
  split <- list(train = train.tb, test = test.tb)
  return(split)
}

assemble.split <- function(split.tb, i) {
  train.nodes <- split.tb$train %>%
    dplyr::filter(fold == i) %>%
    dplyr::pull(nodes) %>%
    sort()
  ntrn <- train.nodes %>%
    length()
  test.nodes <- split.tb$test %>%
    dplyr::filter(fold == i) %>%
    dplyr::pull(nodes) %>%
    sort()
  ntst <- test.nodes %>%
    length()
  return(
    list(
      train = train.nodes,
      test = test.nodes,
      ntrn = ntrn,
      ntst = ntst
    )
  )
}