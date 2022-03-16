assemble.dataset.train <- function(net, train, test, labels, N) {
  train.labs <- labels[train]
  test.labs <- labels[test]
  ec.labs <- c(train.labs, test.labs)
  noec.labs <- c(ec.labs, labels[(N + 1):length(labels)])
  net %>%
    as.data.frame()
  colnames(net) <- labels[1:N]
  rownames(net) <- labels
  net <- net[noec.labs, ec.labs]
  net <- net[noec.labs, train.labs] %>%
    as.matrix()
  return(net)
}

assemble.dataset.test <- function(net, train, test, labels, N) {
  train.labs <- labels[train]
  test.labs <- labels[test]
  ec.labs <- c(train.labs, test.labs)
  noec.labs <- c(ec.labs, labels[(N + 1):length(labels)])
  net %>%
    as.data.frame()
  colnames(net) <- labels[1:N]
  rownames(net) <- labels
  net <- net[noec.labs, ec.labs]
  net <- net[noec.labs, test.labs] %>%
    as.matrix()
  return(net)
}

get.mean.similarity <- function(df, ntrn, ntgt, nsrc=107) {
  net <- df$net
  AIK <- df$aik
  AKI <- df$aki
  AIK[diag(nsrc) == 1] <- NA
  AKI[diag(ntrn) == 1] <- NA
  similarity <- matrix(0, nrow = nsrc, ncol = ntgt)
  # First quadrant ----
  for (i in 1:ntrn) {
    for (j in 1:ntrn) {
      if (i != j) {
        jacp.aik <- AIK[i, j]
        jacp.aki <- AKI[i, j]
        similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
      }
    }
  }
  # Third quadrant ----
  for (i in (ntrn + 1):nsrc) {
    for (j in 1:ntrn) {
      Ni <- which(net[i, ] > 0)
      Nj <- which(net[j, ] > 0)
      NiNj <- Ni[Ni %in% Nj]
      jacp.aik <- AIK[i, j]
      jacp.aki <- mean(AKI[j, NiNj], na.rm = T)
      similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
    }
  }
  # Second quadrant ----
  for (i in 1:ntrn) {
    for (j in (ntrn + 1):ntgt) {
      Ni <- which(net[i, ] > 0)
      Nj <- which(net[j, ] > 0)
      NiNj <- Ni[Ni %in% Nj]
      jacp.aik <- AIK[i, j]
      jacp.aki <- mean(AKI[i, NiNj], na.rm = T)
      similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
    }
  }
  # Forth quadrant ----
  for (i in (ntrn + 1):nsrc) {
    for (j in (ntrn + 1):ntgt) {
      Ni <- which(net[i, ] > 0)
      Nj <- which(net[j, ] > 0)
      NiNj <- Ni[Ni %in% Nj]
      jacp.aik <- AIK[i, j]
      jacp.aki <- mean(AKI[NiNj, NiNj], na.rm = T)
      similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
    }
  }
  return(similarity)
}

col.sum <- function(A, N) {
  x <- c()
  for (i in 1:N) {
    x <- c(x, sum(A[A[, i] > 0, i], na.rm = T))
  }
  return(x)
}

standardize.matrix <- function(matrix) {
  matrix[matrix == 0] <- NA
  mean.matrix <- mean(matrix, na.rm = T)
  sd.matrix <- sd(matrix, na.rm = T)
  matrix <- (matrix -  mean.matrix) / sd.matrix
  matrix[is.na(matrix)] <- 0
  return(matrix)
}

standardize.matrix.reference <- function(matrix, ref) {
  ref[ref == 0] <- NA
  matrix[matrix == 0] <- NA
  mean.ref <- mean(ref, na.rm = T)
  sd.ref <- sd(ref, na.rm = T)
  matrix <- (matrix -  mean.ref) / sd.ref
  matrix[is.na(matrix)] <- 0
  return(matrix)
}

get.dataset <- function(
  net, nt, nodes, labels, split, inst,
  filename="1",save=T) {
  if (save) {
    # Separate training and test set ----
    source("functions/df_to_adj.R")
    ids <- with(net, data.frame(
      source = source, target = target, weight = id)
    ) %>%
      df.to.adj()
    ids <- ids[, 1:nt]
    net <- net %>%
      df.to.adj()
    net <- net[, 1:nt]
    source("functions/get_tracto2016.R")
    distances <- get.tracto2016(labels)
    distances <- distances[, 1:nt]
    distances <- standardize.matrix(distances)
    ## Assembling train data ----
    source("functions/adj_to_df.R")
    id.train <- assemble.dataset.train(
      ids, split$train, split$test, labels, nt
    )
    net.train <- assemble.dataset.train(
      net, split$train, split$test, labels, nt
    )
    dist.train <- assemble.dataset.train(
      distances, split$train, split$test, labels, nt
    ) %>%
      adj.to.df()
    ### Computing AKS ----
    source("functions/compute_aki.R")
    source("functions/compute_aik.R")
    net.aik <- compute.aik(
      net.train %>%
        adj.to.df(),
      nodes
    )
    net.aki <- compute.aki(
      net.train %>%
        adj.to.df(),
      split$ntrn
    )
    ### Get similarities ----
    net.sim <- get.mean.similarity(
      list(
        net = net.train, aik = net.aik, aki = net.aki
      ),
      split$ntrn, nodes
    )
    ### Get similarity train ----
    net.sim.train <- net.sim[, 1:split$ntrn] %>%
      standardize.matrix.reference(net.sim) %>%
        adj.to.df()
    ## Assembling test data ----
    id.test <- assemble.dataset.test(
      ids, split$train, split$test, labels, nt
    )
    net.test <- assemble.dataset.test(
      net, split$train, split$test, labels, nt
    ) %>%
      adj.to.df()
    dist.test <-  assemble.dataset.test(
      distances, split$train, split$test, labels, nt
    )  %>%
      adj.to.df()
    ### Get similarity test and standardize it ----
    net.sim.test <- net.sim[, (split$ntrn+1):nt] %>%
      standardize.matrix.reference(net.sim) %>%
      adj.to.df()
    ## Form train data ----
    net.train <- net.train %>%
      adj.to.df()
    train.data <- data.frame(
      w = net.train$weight,
      sim = net.sim.train$weight,
      dist = dist.train$weight
    )
    matrix.train <- id.train
    id.train <- id.train %>%
      adj.to.df() %>%
      dplyr::filter(weight > 0) %>%
      dplyr::pull(weight)
    ## Form test data ----
    test.data <- data.frame(
      w = net.test$weight,
      sim = net.sim.test$weight,
      dist = dist.test$weight
    )
    matrix.test <- id.test
    id.test <- id.test %>%
      adj.to.df() %>%
      dplyr::filter(weight > 0) %>%
      dplyr::pull(weight)
    ## Save train and test data ----
    df <- list(
      train = train.data,
      test = test.data,
      ids = list(
        train = id.train, test = id.test
      ),
      matids = list(
        train = matrix.train, test = matrix.test
      ),
      split = split
    )
    # saveRDS(
    #   df,
    #   "../RDS/imputation/%s/dataset/%s.rds" %>%
    #     sprintf(inst$common, filename)
    # )
  } else
    df <- readRDS(
      "../RDS/imputation/%s/dataset/%s.rds" %>%
        sprintf(inst$common, filename)
    )
  data <- structure(
    list(
      data = df$train %>%
        rbind(df$test) %>%
        dplyr::as_tibble(),
      in_id = nrow(df$train) %>%
        seq_len(),
      out_id = (nrow(df$train) + 1):(nrow(df$train) + nrow(df$test))
    ),
    class = "rsplit"
  )
  return(list(
    data = data, ids = df$ids, matids = df$matids)
  )
}