#' Relative mean absolute error
#'
#' Calculate the mean absolute error. This metric is in the same units as the
#' original data.
#'
#' @family numeric metrics
#' @family accuracy metrics
#' @templateVar metric_fn rmae
#' @template return
#'
#' @inheritParams rmse
#'
#' @author Jorge Martinez using code from Max Kuhn
#'
#' @template examples-numeric
#'
#' @export

rmae <- function(data, ...) {
  UseMethod("rmae")
}
rmae <- yardstick::new_numeric_metric(
  rmae,
  direction = "minimize"
)

#' @rdname rmae
#' @export
rmae.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  yardstick::metric_summarizer(
    metric_nm = "rmae",
    metric_fn = rmae_vec,
    data = data,
    truth = !!dplyr::enquo(truth),
    estimate = !!dplyr::enquo(estimate),
    na_rm = na_rm
  )
  
}

#' @export
#' @rdname rmae
rmae_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  
  rmae_impl <- function(truth, estimate) {
    mean( abs((truth - estimate)/truth))
  }
  
  yardstick::metric_vec_template(
    metric_impl = rmae_impl,
    truth = truth,
    estimate = estimate,
    na_rm = na_rm,
    cls = "numeric"
  )
  
}