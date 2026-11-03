#' Low-Variability Periods
#'
#' Extracts the periods (hours of the day) of low variability of a specific parameter. Optionally crops the dataset to these periods.
#'
#'
#' @details
#' Let \eqn{X_{it}} denote the measurements at a location \eqn{i} at each timestamp \eqn{t \in \mathcal{T}}, and \eqn{h_t \in \{0, 1, \dots, 23\}} the hour of the day of \eqn{t}.
#'
#' All measurements recorded at hour \eqn{h} are grouped as:
#' \deqn{
#' X_h = \{ X_{it} \mid h_t = h, \, \forall i, t \}.
#' }{X_h = all X_it with hour h}
#'
#' Missing values are excluded from the calculation. The variance of the measurements at each hour \eqn{h} is calculated as:
#' \deqn{
#' \sigma_h^2 = \frac{1}{N_h - 1} \sum_{X \in X_h} (X - \overline{X}_h)^2,
#' }{var_h = sum((X - mean_h)^2) / (N_h - 1)}
#' \itemize{
#'   \item \eqn{N_h = |X_h|}: number of available (non-missing) measurements at hour \eqn{h}.
#'   \item \eqn{\overline{X}_h = \frac{1}{N_h} \sum_{X \in X_h} X}: mean value at hour \eqn{h}.
#' }
#'
#' The set of the \eqn{n_{\text{min}}}{n_min} (`nmin`) hours with the smallest variances is defined as:
#' \deqn{
#' H^{\ast} = \{ h_1, h_2, \dots, h_{n_{\text{min}}} \}, \quad \sigma_{h_1}^2 \leq \sigma_{h_2}^2 \leq \dots \leq \sigma_{h_{n_{\text{min}}}}^2,
#' }{H* = the n_min hours with the smallest variance}
#' where \eqn{h_1} is the hour with the lowest variance, \eqn{h_2} the second lowest, and so on.
#'
#' The timestamps of the low-variability periods are then:
#' \deqn{
#' \mathcal{T}^{\ast} = \{ t \in \mathcal{T} \mid h_t \in H^{\ast} \}.
#' }{T* = timestamps t with hour in H*}
#'
#' @param data A data frame containing the dataset to be used.
#' @param date.column character. Name of the date/time column in `data`.
#' @param data.column character. Name of the column in `data` used to calculate the variance.
#' @param nmin numeric or `NULL`. Number of hours with the minimum variability. If `NULL`, it is set to `24`. Default is `NULL`.
#' @param crop logical. If `TRUE`, `data` is cropped to include the periods with low variability. Default is `TRUE`.
#
#' @return Returns a data frame (subset of `data`) including only the periods with low variability when `crop = TRUE`. If `crop = FALSE`, it returns a numeric vector representing the low-variability periods.
#'
#' @seealso \code{\link{correction_process}}
#'
#' @export
min_variability_function <- function(data, date.column, data.column, nmin, crop = TRUE) {
  if (is.null(nmin)) {
    nmin <- 24
  }
  stopifnot(is.numeric(nmin), length(nmin) == 1, !is.na(nmin), nmin >= 1, nmin <= 24)
  #
  data.diurnal <- data %>%
    mutate(Hour = hour(.data[[date.column]])) %>%
    group_by(.data[['Hour']]) %>%
    summarise(variability = var(.data[[data.column]], na.rm = TRUE), .groups = "drop") %>%
    arrange(.data[['variability']])
  min.var <- sort(slice_head(data.diurnal, n = nmin)[['Hour']])
  if (!crop) {
    return(min.var)
  } else {
    data <- data %>%
      mutate(Hour = hour(.data[[date.column]])) %>%
      filter(.data[['Hour']] %in% min.var) %>%
      select(-all_of("Hour"))
    return(data)
  }
}
