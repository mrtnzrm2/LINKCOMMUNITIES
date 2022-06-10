rmae <- function(truth, pred) {
  if (0 %in% truth) {
    pred <- pred[truth != 0]
    truth <- truth[truth != 0]
  }
  return(mean(abs((truth - pred) / truth)))
}