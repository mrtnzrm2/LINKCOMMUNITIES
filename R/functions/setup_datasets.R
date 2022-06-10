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
  similarity <- matrix(0, nrow = nsrc, ncol = ntgt)
  # First quadrant ----
  mean_exact_sim <- array(1, dim = c(ntrn, ntrn, 2))
  mean_exact_sim[, , 1] <- AIK[1:ntrn, 1:ntrn]
  mean_exact_sim[, , 2] <- AKI[1:ntrn, 1:ntrn]
  similarity[1:ntrn, 1:ntrn] <- mean_exact_sim %>%
    apply(c(1, 2), mean)
  # Third quadrant ----
  for (i in (ntrn + 1):nsrc) {
    for (j in 1:ntrn) {
      Ni <- which(net[i, ] > 0)
      Nj <- which(net[j, ] > 0)
      NiNj <- Ni[Ni %in% Nj]
      if (length(NiNj) == 0)
          next
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
      if (length(NiNj) == 0)
          next
      jacp.aik <- AIK[i, j]
      jacp.aki <- mean(AKI[i, NiNj], na.rm = T)
      similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
    }
  }
  # Forth quadrant ----
  for (i in (ntrn + 1):nsrc) {
    for (j in (ntrn + 1):ntgt) {
      if (i != j) {
        Ni <- which(net[i, ] > 0)
        Nj <- which(net[j, ] > 0)
        NiNj <- Ni[Ni %in% Nj]
        if (length(NiNj) == 0)
          next
        jacp.aik <- AIK[i, j]
        jacp.aki <- mean(AKI[NiNj, NiNj], na.rm = T)
        similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
      }
    }
  }
  return(similarity)
}

get.mean.similarity_2 <- function(df, ntrn, ntgt, nsrc=107) {
  net <- df$net
  AIK <- df$aik
  AKI <- df$aki
  similarity <- matrix(0, nrow = nsrc, ncol = ntgt)
  # First quadrant ----
  mean_exact_sim <- array(1, dim = c(ntrn, ntrn, 2))
  mean_exact_sim[, , 1] <- AIK[1:ntrn, 1:ntrn]
  mean_exact_sim[, , 2] <- AKI[1:ntrn, 1:ntrn]
  similarity[1:ntrn, 1:ntrn] <- mean_exact_sim %>%
    apply(c(1, 2), mean)
  # Third quadrant ----
  for (i in (ntrn + 1):nsrc) {
    for (j in 1:ntrn) {
      Ni <- which(net[i, ] > 0)
      Nj <- which(net[j, ] > 0)
      NiNj <- Ni[Ni %in% Nj]
      if (length(NiNj) == 0) {
        similarity[i, j] <- 0.043
        next
      }
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
      if (length(NiNj) == 0) {
        similarity[i, j] <- 0.043
        next
      }
      jacp.aik <- AIK[i, j]
      jacp.aki <- mean(AKI[i, NiNj], na.rm = T)
      similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
    }
  }
  # Forth quadrant ----
  for (i in (ntrn + 1):nsrc) {
    for (j in (ntrn + 1):ntgt) {
      if (i != j) {
        Ni <- which(net[i, ] > 0)
        Nj <- which(net[j, ] > 0)
        NiNj <- Ni[Ni %in% Nj]
        if (length(NiNj) == 0) {
          similarity[i, j] <- 0.023
          next
        }
        jacp.aik <- AIK[i, j]
        jacp.aki <- mean(AKI[NiNj, NiNj], na.rm = T)
        similarity[i, j] <- mean(c(jacp.aki, jacp.aik))
      }
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

get.dataset.fln_dist <- function(
  net, nt, nodes, labels, split, inst) {
  # Separate training and test set ----
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/standardized.R")
  #*** Ids from nodesxnt matrix v
  ids <- with(net, data.frame(
    source = source, target = target, weight = id)
  ) %>%
    #*** to nodesxnt matrix v
    df.to.adj()
  ids <- ids[, 1:nt]
  #*** net here is a data frame v
  net <- net %>%
    #*** to nodesxnt matrix v
    df.to.adj()
  #*** nodesxnt matrix v
  net <- net[, 1:nt]
  source("functions/get_tracto2016.R")
  distances <- get.tracto2016(labels)
  #*** self-loops na introduced v
  distances <- standardized.diag(distances, nodes)
  #*** nodesxnt size matrix v
  distances <- distances[, 1:nt]
  ## Assembling train data ----
  #*** nodesxntrn size matrix v
  id_train <- assemble.dataset.train(
    ids, split$train, split$test, labels, nt
  )
  #*** nodesxntrn size matrix v
  dist_train <- assemble.dataset.train(
    distances, split$train, split$test, labels, nt
  ) %>%
    #*** to data frame v
    adj.to.df()
  ## Assembling test data ----
  #*** nodesx(nt-ntrn) size matrix
  id_test <- assemble.dataset.test(
    ids, split$train, split$test, labels, nt
  )
  #*** nodesx(nt-ntrn) size matrix v
  dist_test <-  assemble.dataset.test(
    distances, split$train, split$test, labels, nt
  )  %>%
    #*** to data frame v
    adj.to.df()
  # Get train and test set data and features ----
  assembled_data <- assemble.dataset.fln_dist(
    net, split$train, split$test, labels, nt
  )
  #*** train and test weights as list v
  #*** zeros and self-loops included v
  train_weight <- assembled_data$train
  test_weight <- assembled_data$test
  #*** Features matrices v
  features_train <- assembled_data$features_train
  features_test <- assembled_data$features_test
  # Build up train and test datasets ----
  #*** train set v
  data.train <- train_weight %>%
    dplyr::bind_cols(dist_train$weight)
  colnames(data.train) <- c("w", "dist")
  data.train <- data.train %>%
    dplyr::bind_cols(features_train)
  #*** Take out all zeros from train set v
  data.train <- data.train[data.train$w > 0, ]
  # Train ids as matrix and list ----
  matrix_train <- id_train
  id_train <- id_train %>%
    adj.to.df() %>%
    dplyr::filter(weight > 0) %>%
    dplyr::pull(weight)
  #*** Test set v
  data.test <- test_weight %>%
    dplyr::bind_cols(dist_test$weight)
  colnames(data.test) <- c("w", "dist")
  data.test <- data.test %>%
    dplyr::bind_cols(features_test)
  #*** Take out self_loops from test set v
  self_loops <- is.na(dist_test$weight)
  data.test <- data.test[!self_loops, ]
  # Test ids as matrix and list ----
  matrix_test <- id_test
  id_test <- id_test %>%
    adj.to.df() %>%
    dplyr::pull(weight)
  id_test <- id_test[!self_loops]
  ## Save train and test data ----
  data <- structure(
    list(
      data = data.train %>%
        rbind(data.test) %>%
        dplyr::as_tibble(),
      in_id = nrow(data.train) %>%
        seq_len(),
      out_id = (nrow(data.train) + 1):(nrow(data.train) + nrow(data.test))
    ),
    class = "rsplit"
  )
  return(list(
      data = data,
      ids = list(
        train = id_train, test = id_test
      ),
      matids = list(
        train = matrix_train, test = matrix_test
      )
    )
  )
}

get.dataset.sim_dist <- function(
  net, nt, nodes, labels, split, inst, zero=F) {
  ## Train set does not have any zero
  ## Test set does not have self-loops
  # Separate training and test set ----
  source("functions/df_to_adj.R")
  source("functions/standardized.R")
  #*** Ids from nodesxnt matrix v
  ids <- with(net, data.frame(
    source = source, target = target, weight = id)
  ) %>%
    #*** to nodesxnt matrix v
    df.to.adj()
  ids <- ids[, 1:nt]
  #*** net here is a data frame v
  net <- net %>%
    #*** to nodesxnt matrix v
    df.to.adj()
  #*** nodesxnt matrix v
  net <- net[, 1:nt]
  source("functions/get_tracto2016.R")
  distances <- get.tracto2016(labels)
  #*** self-loops na introduced v
  distances <- divide_max(distances, nodes)
  #*** nodesxnt size matrix v
  distances <- distances[, 1:nt]
  ## Assembling train data ----
  source("functions/adj_to_df.R")
  #*** nodesxntrn size matrix v
  id_train <- assemble.dataset.train(
    ids, split$train, split$test, labels, nt
  )
  #*** nodesxntrn size matrix v
  net_train <- assemble.dataset.train(
    net, split$train, split$test, labels, nt
  )
  #*** nodesxntrn size matrix v
  dist_train <- assemble.dataset.train(
    distances, split$train, split$test, labels, nt
  ) %>%
    #*** to data frame v
    adj.to.df()
  ### Computing AKS ----
  source("functions/compute_aki.R")
  source("functions/compute_aik.R")
  #*** nodesxnodes size matrix v
  net_aik <- compute.aik(
    net_train %>%
      adj.to.df(),
    nodes
  )
  #*** ntrnxntrn size matrix v
  net_aki <- compute.aki(
    net_train %>%
      adj.to.df(),
    split$ntrn
  )
  ### Get similarities ----
  #*** nodesxnt size matrix v
  net_sim <- get.mean.similarity_2(
    list(
      net = net_train, aik = net_aik, aki = net_aki
    ),
    split$ntrn, nt
  ) 
  # %>%
  #   #***  block diag self-loops na introduced v
  #   standardized.block.diag()
  ### Get similarity train ----
  #*** nodesxntrn size matrix v
  net_sim_train <- net_sim[, 1:split$ntrn] %>%
      adj.to.df()
  ## Assembling test data ----
  #*** nodesx(nt-ntrn) size matrix
  id_test <- assemble.dataset.test(
    ids, split$train, split$test, labels, nt
  )
  #*** nodesx(nt-ntrn) size matrix v
  net_test <- assemble.dataset.test(
    net, split$train, split$test, labels, nt
  ) %>%
    #*** to data frame v
    adj.to.df()
  #*** nodesx(nt-ntrn) size matrix v
  dist_test <-  assemble.dataset.test(
    distances, split$train, split$test, labels, nt
  )  %>%
    #*** to data frame v
    adj.to.df()
  ### Get similarity test and standardize it ----
  #*** nodesx(nt-ntrn) size matrix v
  net_sim_test <- net_sim[, (split$ntrn + 1):nt] %>%
    #*** to data frame v
    adj.to.df()
  ## Form train data ----
  #*** nodesxntrn size matrix v
  net_train <- net_train %>%
    #*** to data frame v
    adj.to.df()
  # With zeros or without zeros ----
  obj <- get_zeros(
    net_train, net_sim_train, dist_train, id_train,
    net_test, net_sim_test, dist_test, id_test,
    split, on = zero
  )
  data <- obj$data
  train_data <- data$train
  test_data <- data$test
  id <- obj$id
  id_train <- id$train
  id_test <- id$test
  mat <- obj$matrix
  matrix_train <- mat$train
  matrix_test <- mat$test
  st <- obj$st
  train_st <- st$train
  test_st <- st$test
  ## Save train and test data ----
  data <- structure(
    list(
      data = train_data %>%
        rbind(test_data) %>%
        dplyr::as_tibble(),
      in_id = nrow(train_data) %>%
        seq_len(),
      out_id = (nrow(train_data) + 1):(nrow(train_data) + nrow(test_data))
    ),
    class = "rsplit"
  )
  return(list(
      data = data,
      ids = list(
        train = id_train, test = id_test
      ),
      matids = list(
        train = matrix_train, test = matrix_test
      ),
      st = list(
        train = train_st,
        test = test_st
      ),
      sim_df = list(
        train = net_sim_train %>%
          df.to.adj(),
        test = net_sim_test %>%
          df.to.adj()
      )
    )
  )
}

get_zeros <- function(
  net_train, net_sim_train, dist_train, id_train,
  net_test, net_sim_test, dist_test, id_test,
  split, on=T
  ) {
  if (on) {
    train_data <- data.frame(
      w = net_train$weight,
      sim = net_sim_train$weight,
      dist = dist_train$weight
    )
    train_st <- data.frame(
      source = net_sim_train$source,
      target = net_sim_train$target
    )
    # Train ids as matrix and list ----
    matrix_train <- id_train
    id_train <- id_train %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    ## Form test data ----
    test_data <- data.frame(
      w = net_test$weight,
      sim = net_sim_test$weight,
      dist = dist_test$weight
    )
    test_st <- data.frame(
      source = net_sim_test$source,
      target = net_sim_test$target
    )
    # Test ids as matrix and list ----
    matrix_test <- id_test
    id_test <- id_test %>%
      adj.to.df() %>%
      dplyr::pull(weight)
  } else {
    # Take out zeros ----
    zeros <- net_train$w == 0
    train_data <- data.frame(
      w = net_train$weight[!zeros],
      sim = net_sim_train$weight[!zeros],
      dist = dist_train$weight[!zeros]
    )
    train_st <- data.frame(
      source = net_sim_train$source[!zeros],
      target = net_sim_train$target[!zeros]
    )
    # Train ids as matrix and list ----
    matrix_train <- id_train
    id_train <- id_train %>%
      adj.to.df() %>%
      dplyr::filter(weight > 0) %>%
      dplyr::pull(weight)
    ## Form test data ----
    ### Take out self-loops from test set ----
    self_loops <- is.na(net_sim_test$weight)
    test_data <- data.frame(
      w = net_test$weight[!self_loops],
      sim = net_sim_test$weight[!self_loops],
      dist = dist_test$weight[!self_loops]
    )
    test_st <- data.frame(
      source = net_sim_test$source[!self_loops],
      target = net_sim_test$target[!self_loops] + split$ntrn
    )
    # Test ids as matrix and list ----
    matrix_test <- id_test
    id_test <- id_test %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    id_test <- id_test[!self_loops]
  }
  return(
    list(
      data = list(
        train = train_data,
        test = test_data
      ),
      id = list(
        train = id_train,
        test = id_test
      ),
      matrix = list(
        train = matrix_train,
        test = matrix_test
      ),
      st = list(
        train = train_st,
        test = test_st
      )
    )
  )
}

assemble.dataset.fln_dist <- function(net, train, test, labels, N) {
  train_labs <- labels[train]
  test_labs <- labels[test]
  ec_labs <- c(train_labs, test_labs)
  noec_labs <- c(ec_labs, labels[(N + 1):length(labels)])
  net %>%
    as.data.frame()
  colnames(net) <- labels[1:N]
  rownames(net) <- labels
  net <- net[noec_labs, ec_labs]
  net_train <- net[, train_labs] %>%
    adj.to.df()
  net_test <- net[, test_labs] %>%
    adj.to.df()
  ntrn <- train %>%
    length()
  train_features <- matrix(
    0, nrow(net_train), 2 * ntrn
  )
  test_features <- matrix(
    0, nrow(net_test), 2 * ntrn
  )
  for (i in seq_len(nrow(net_train))) {
    train_features[i, 1:ntrn] <- net[
      net_train$source[i], seq_len(ntrn)
    ]
    train_features[i, (ntrn + 1):(2 * ntrn)] <- net[
      net_train$target[i], seq_len(ntrn)
    ]
  }
  for (i in seq_len(nrow(net_test))) {
    test_features[i, 1:ntrn] <- net[
      net_test$source[i], seq_len(ntrn)
    ]
    test_features[i, (ntrn + 1):(2 * ntrn)] <- net[
      net_test$target[i] + ntrn, seq_len(ntrn)
    ]
  }
  return(
    list(
      train = net_train %>%
        dplyr::pull(weight),
      test = net_test %>%
        dplyr::pull(weight),
      features_train = train_features %>%
        dplyr::as_tibble(),
      features_test = test_features %>%
        dplyr::as_tibble()
    )
  )
}