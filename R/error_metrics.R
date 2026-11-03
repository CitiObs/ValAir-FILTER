#' Error Metrics
#'
#' Computes common statistical error indicators between two numeric vectors: the mean bias difference (MBD), the mean absolute difference (MAD), the root mean square difference (RMSD), the mean absolute percentage error (MAPE) and the normalized mean standard deviation (nMSD).
#'
#'
#' @details
#' Let \eqn{X_i} denote the predicted values and \eqn{Y_i} the observed values, with \eqn{i = 1, \dots, n}. Missing values are removed before the calculation, unless `all = TRUE`. The following metrics are included:
#'
#' \strong{Mean Bias Difference (MBD) and Mean Absolute Difference (MAD)}
#'
#' The MBD is calculated by `mbd()` as:
#' \deqn{
#' MBD = \frac{1}{n} \sum_{i=1}^n (X_{i} - Y_{i}),
#' }{MBD = mean(X - Y)}
#' and, if `absolute = TRUE`, the MAD as:
#' \deqn{
#' MAD = \frac{1}{n} \sum_{i=1}^n |X_{i} - Y_{i}|.
#' }{MAD = mean(|X - Y|)}
#'
#' \strong{Root Mean Square Difference (RMSD)}
#'
#' The RMSD is calculated by `rmsd()` as:
#' \deqn{
#' RMSD = \sqrt{ \frac{1}{n} \sum_{i=1}^n (X_{i} - Y_{i})^2 }.
#' }{RMSD = sqrt(mean((X - Y)^2))}
#'
#' \strong{Normalization of MBD and RMSD}
#'
#' If `norm = TRUE`, the normalized metrics are calculated by dividing by a measure of the observed values, according to `scale`:
#' \deqn{
#' nMBD = \frac{MBD}{\bar{Y}}, \quad nRMSD = \frac{RMSD}{\bar{Y}} \quad (\text{`scale = FALSE`}),
#' }{nMBD = MBD / mean(Y), nRMSD = RMSD / mean(Y)}
#' \deqn{
#' nMBD = \frac{MBD}{sd(Y)}, \quad nRMSD = \frac{RMSD}{sd(Y)} \quad (\text{`scale = TRUE`}),
#' }{nMBD = MBD / sd(Y), nRMSD = RMSD / sd(Y)}
#' where \eqn{\bar{Y}} is the mean and \eqn{sd(Y)} the standard deviation of the observed values.
#'
#' \strong{Mean Absolute Percentage Error (MAPE)}
#'
#' The MAPE is calculated by `mape()` as:
#' \deqn{
#' MAPE = \frac{1}{n} \sum_{i=1}^n \left| \frac{X_{i} - Y_{i}}{Y_i} \right|.
#' }{MAPE = mean(|(X - Y) / Y|)}
#'
#' \strong{Normalized Mean Standard Deviation (nMSD)}
#'
#' The nMSD is calculated by `nmsd()` as:
#' \deqn{
#' nMSD = \frac{sd(X) - sd(Y)}{sd(Y)}.
#' }{nMSD = (sd(X) - sd(Y)) / sd(Y)}
#'
#' These metrics are useful for evaluating the accuracy, the bias and the variability of the predictions compared to the observations.
#'
#' @param x numeric. Vector of predicted values.
#' @param y numeric. Vector of observed values.
#' @param all logical. If `FALSE`, missing values are removed before the calculation. Default is `FALSE`.
#' @param norm logical. If `TRUE`, the normalized (relative) metric is calculated. Applicable only to `mbd()` and `rmsd()`. Default is `FALSE`.
#' @param scale logical. If `FALSE`, the metric is normalized by the mean of the observed values; if `TRUE`, by their standard deviation. Used only if `norm = TRUE`. Default is `FALSE`.
#' @param absolute logical. If `TRUE`, the absolute differences are used (MAD instead of MBD). Applicable only to `mbd()`. Default is `FALSE`.
#'
#' @name error_metrics
#' @aliases mbd rmsd mape nmsd
#' @export mbd
#' @export rmsd
#' @export mape
#' @export nmsd
mbd <- function(x, y, all = FALSE, norm = FALSE, scale = FALSE, absolute = FALSE) {
  xy <- data.frame(x, y)
  if (!all) xy <- xy %>% drop_na()
  x <- xy[['x']]
  y <- xy[['y']]
  mb <- if (!absolute) mean(x - y, na.rm = TRUE) else mean(abs(x - y), na.rm = TRUE)
  if (norm) {
    mb <- if (!scale) mb/mean(y, na.rm = TRUE) else mb/sd(y, na.rm = TRUE)
  }
  return(mb)
}

#' @rdname error_metrics
#' @export
rmsd <- function(x, y, all = FALSE, norm = FALSE, scale = FALSE) {
  xy <- data.frame(x, y)
  if (!all) xy <- xy %>% drop_na()
  x <- xy[['x']]
  y <- xy[['y']]
  rms <- sqrt(mean((x - y)^2, na.rm = TRUE))
  if (norm) {
    rms <- if (!scale) rms/mean(y, na.rm = TRUE) else rms/sd(y, na.rm = TRUE)
  }
  return(rms)
}

#' @rdname error_metrics
#' @export
mape <- function(x, y, all = FALSE) {
  xy <- data.frame(x, y)
  if (!all) xy <- xy %>% drop_na()
  x <- xy[['x']]
  y <- xy[['y']]
  mpe <- mean(abs(x - y) / abs(y), na.rm = TRUE)
  return(mpe)
}

#' @rdname error_metrics
#' @export
nmsd <- function(x, y, all = FALSE) {
  xy <- data.frame(x, y)
  if (!all) xy <- xy %>% drop_na()
  x <- xy[['x']]
  y <- xy[['y']]
  nms <- (sd(x, na.rm = TRUE) - sd(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
  return(nms)
}
