#' Temporal Outlier Detection Test
#'
#' Detects abnormal values (outliers) within a temporal window, and assigns a quality flag. Optionally creates a column of valid values.
#'
#'
#' @details
#' Let \eqn{X_t} denote the measurements at a location. For each time step
#' \eqn{t_i}, the temporal window \eqn{\mathcal{T}} is defined as:
#' \deqn{
#' \mathcal{T} = [t_i - D, t_i],
#' }
#' where \eqn{D} (`outlier.duration`) is the duration of the window.
#'
#' Missing values are excluded from the calculation. Let \eqn{\mathcal{T}^{\ast} \subseteq \mathcal{T}} be the set of available
#' (non-missing) values in the window, and \eqn{n_{\mathcal{T}^{\ast}}} their number. The statistical index \eqn{Z_t} for detecting outliers is calculated as:
#' \deqn{
#' Z_t = \frac{|X_t - m_{\mathcal{T}^{\ast}}|}{\text{MAD}_{\mathcal{T}^{\ast}}},
#' }
#' \itemize{
#'   \item \eqn{m_{\mathcal{T}^{\ast}}}: median of the available (non-missing) values within \eqn{\mathcal{T}^{\ast}}.
#'   \item \eqn{\text{MAD}_{\mathcal{T}^{\ast}}}: median absolute deviation, defined as:
#'   \deqn{
#'   \text{MAD}_{\mathcal{T}^{\ast}} = c \cdot \text{median}_{t \in \mathcal{T}^{\ast}} \left( |X_t - m_{\mathcal{T}^{\ast}}| \right)
#'   }
#'   \item \eqn{c}: scale factor (for normally distributed data, \eqn{c = 1.4826}; here \eqn{c = 1}).
#' }
#'
#' The quality flag \eqn{f_{t}^{\text{outlier}}} is then defined as:
#' \deqn{
#' f_{t}^{\text{outlier}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: } n_{\mathcal{T}^{\ast}} < n_{min}, \\
#' 1 & \text{Outlier: } n_{\mathcal{T}^{\ast}} \geq n_{min} \text{ and } Z_t \geq q, \\
#' 2 & \text{Non-outlier: } n_{\mathcal{T}^{\ast}} \geq n_{min} \text{ and } Z_t < q,
#' \end{cases}
#' }
#' \itemize{
#'   \item \eqn{n_{min}}{n_min} (`compl.outlier.duration`): minimum number of non-missing values required to evaluate the test, given either as a number of measurements or as a fraction \eqn{p} of the window duration \eqn{D} (`outlier.duration`), in which case \eqn{n_{min} = p D}{n_min = p * D}.
#'   \item \eqn{q} (`outlier.threshold`): threshold for the statistical index.
#' }
#'
#' @param data A data frame containing the dataset to be tested
#' @param date.column character. Name of the date/time column in `data`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location ID. If `NULL`, the test is applied to the whole dataset; otherwise, it is applied to each individual location. Default is `NULL`.
#' @param outlier.duration numeric. Length of the rolling window.
#' @param compl.outlier.duration numeric. Minimum data completeness required within `outlier.duration` to run the test, expressed either as a number of measurements or as a fraction of `outlier.duration`.
#' @param outlier.threshold numeric. Threshold of the calculated statistical index. Default is `10`.
#' @param new logical. If TRUE, adds a column `<param>_outlier` containing the valid values; flagged values are replaced with `NA`. Default is `FALSE`.
#' @param param character. Base name of the new column (e.g. `PM2.5` gives `PM2.5_outlier`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_outlier}: integer. Quality flags based on the statistical test, \code{0} = insufficient data, \code{1} = outlier, \code{2} = non-outlier. It is `NA` where the tested value is missing.
#'   \item \code{<param>_outlier}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#'
#' @export
outlier_qctest <- function(data,
                           date.column,
                           data.column,
                           id.column = NULL,
                           outlier.duration,
                           compl.outlier.duration,
                           outlier.threshold,
                           new = FALSE,
                           param = 'PM2.5') {
  # Argument checks
  stopifnot(is.character(date.column), length(date.column) == 1, date.column %in% colnames(data))
  stopifnot(is.character(data.column), length(data.column) == 1, data.column %in% colnames(data))
  stopifnot(outlier.duration > 0, !is.null(outlier.duration))
  stopifnot(compl.outlier.duration >= 0, !is.null(compl.outlier.duration))
  stopifnot(outlier.threshold >= 0)

    stopifnot(!missing(param) && !is.null(param))
  if (!is.null(id.column)) {
    stopifnot(is.character(id.column), length(id.column) == 1, id.column %in% colnames(data))
  }
  #
  compl.outlier.duration <- ifelse(compl.outlier.duration < 1,
                                   round(compl.outlier.duration * outlier.duration),
                                   compl.outlier.duration)
  #
  outlier_fun <- function(df) {
    ind.nona <- which(!is.na(df[[data.column]]))
    if (length(ind.nona) == 0) {
      df <- df %>% mutate(flag_outlier = NA_real_)
    } else {
      start.nona <- df[[date.column]][ind.nona[1]]
      end.nona <- df[[date.column]][ind.nona[length(ind.nona)]]
      df.sub <- df %>%
        filter(.data[[date.column]] >= start.nona, .data[[date.column]] <= end.nona)
      #
      roll.diff.med <- zoo::rollapply(df.sub[[data.column]], outlier.duration, FUN = function(x) {
        med <- if (all(is.na(x))) NA_real_ else round(median(x, na.rm = TRUE), 1)
        if (is.na(med) || is.infinite(med)) NA_real_ else abs(round(x[length(x)] - med, 1))
      }, by.column = FALSE, fill = NA_real_, align = 'right')
      roll.mad <- zoo::rollapply(df.sub[[data.column]], outlier.duration, FUN = function(x) {
        mad.x <- if (all(is.na(x))) NA_real_ else round(mad(x, na.rm = TRUE, constant = 1), 1)
        if (is.na(mad.x) || is.infinite(mad.x)) NA else mad.x
      }, by.column = FALSE, fill = NA_real_, align = 'right')
      roll.nobs <- zoo::rollapply(df.sub[[data.column]], outlier.duration, FUN = function(x) {
        sum(!is.na(x))
      }, by.column = FALSE, fill = NA_real_, align = 'right')
      #
      df.sub <- df.sub %>%
        mutate(flag_outlier = case_when(
          roll.nobs < compl.outlier.duration ~ 0,
          roll.diff.med >= outlier.threshold * roll.mad & roll.nobs >= compl.outlier.duration ~ 1,
          roll.diff.med < outlier.threshold * roll.mad & roll.nobs >= compl.outlier.duration ~ 2,
          TRUE ~ NA_real_
        ))
      df.out <- df %>%
        filter(.data[[date.column]] < start.nona | .data[[date.column]] > end.nona)
      df <- bind_rows(df.out, df.sub) %>%
        arrange(.data[[date.column]]) %>%
        mutate(flag_outlier = ifelse(is.na(.data[[data.column]]), NA_real_, .data[["flag_outlier"]]))
    }
    #
    if (new) {
      df <- df %>%
        mutate(!!paste0(param, "_outlier") := if_else(
          .data[["flag_outlier"]] %in% c(0, 1) | is.na(.data[["flag_outlier"]]),
          NA_real_,
          .data[[data.column]]
        ))
    }
    return(df)
  }
  #
  if (!is.null(id.column)) {
    data <- data %>%
      group_by(across(all_of(id.column))) %>%
      group_map(~outlier_fun(.x), .keep = TRUE) %>%
      bind_rows()
  } else {
    data <- outlier_fun(data)
  }
  return(as_tibble(data))
}
