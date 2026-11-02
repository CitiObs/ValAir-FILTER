#' Constant Value Test
#'
#' Detects flatline behavior, i.e. consecutive measurements that remain almost constant within a temporal window, and assigns a quality flag. Optionally creates a column of valid values.
#'
#'
#' @details
#' Let \eqn{X_t} denote the measurements at a location. For each time step
#' \eqn{t_i}, the temporal window \eqn{\mathcal{T}} is defined as:
#' \deqn{
#' \mathcal{T} = [t_i - D, t_i],
#' }
#' where \eqn{D} (`persist.duration`) is the duration of the window.
#'
#' Missing values are excluded from the calculation. Let \eqn{\mathcal{T}^{\ast} \subseteq \mathcal{T}} be the set of available
#' (non-missing) values in the window, and \eqn{n_{\mathcal{T}^{\ast}}} their number. The variability \eqn{\text{S}_{\mathcal{T}^{\ast}}} within the window is measured by:
#
#' \itemize{
#'   \item Range (`method = range`):
#'   \deqn{
#'   \text{Range}_{\mathcal{T}^{\ast}} = \max_{t \in \mathcal{T}^{\ast}}(X_t) - \min_{t \in \mathcal{T}^{\ast}}(X_t)
#'   }
#'   \item Standard Deviation (`method = range`):
#'   \deqn{
#'   \text{SD}_{\mathcal{T}^{\ast}} = \sqrt{\frac{1}{n_{\mathcal{T}^{\ast}} - 1} \sum_{t \in \mathcal{T}^{\ast}} (X_t - E_{\mathcal{T}^{\ast}})^2}
#'   }
#'   \deqn{
#'   E_{\mathcal{T}^{\ast}} = \frac{1}{n_{\mathcal{T}^{\ast}}} \sum_{t \in \mathcal{T}^{\ast}} X_t
#'   }
#' }
#'
#' The quality flag \eqn{f_{t}^{\text{constant}}} is then defined as:
#' \deqn{
#' f_{t}^{\text{constant}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: } n_{\mathcal{T}^{\ast}} < n_{min}, \\
#' 1 & \text{Constant: } n_{\mathcal{T}^{\ast}} \geq n_{min} \text{ and } \text{S}_{\mathcal{T}^{\ast}} \leq \text{S}_{min}, \\
#' 2 & \text{Non-constant: } n_{\mathcal{T}^{\ast}} \geq n_{min} \text{ and } \text{S}_{\mathcal{T}^{\ast}} > \text{S}_{min},
#' \end{cases}
#' }
#' \itemize{
#'   \item \eqn{n_{min}}{n_min} (`compl.persist.duration`): minimum number of non-missing values required to evaluate the test, given either as a number of measurements or as a fraction \eqn{p} of the window duration \eqn{D} (`persist.duration`), in which case \eqn{n_{min} = p D}{n_min = p * D}.
#'   \item \eqn{\mathrm{S}_{min}}{S_min} (`min.variation`): minimum acceptable threshold for the statistical index within the window.
#' }
#'
#' @param data A data frame containing the dataset to be tested
#' @param date.column character. Name of the date/time column in `data`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location ID. If `NULL`, the test is applied to the whole dataset; otherwise, it is applied to each individual location. Default is `NULL`.
#' @param persist.duration numeric. Length of the rolling window.
#' @param compl.persist.duration numeric. Minimum data completeness required within `persist.duration` to run the test, expressed either as a number of measurements or as a fraction of `persist.duration`.
#' @param method character. Statistical index used to calculate the variability. Either `range` or `std`.
#' @param min.variation numeric. Minimum range below which values are considered constant. Default is `0.1`.
#' @param new logical. If `TRUE`, adds a column `<param>_constant_value` containing the valid values; flagged values are replaced with with `NA`. Default is `FALSE`.
#' @param param character. Base name of the new column (e.g. `PM2.5` gives `PM2.5_constant_value`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_constant_value}: integer. Quality flags based on the statistical test, \code{0} = insufficient data, \code{1} = constant, \code{2} = non-constant. It is `NA` where the tested value is missing.
#'   \item \code{<param>_constant_value}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
constant_value_qctest <- function(data,
                                  date.column,
                                  data.column,
                                  id.column = NULL,
                                  persist.duration,
                                  compl.persist.duration,
                                  method = c('std','range'),
                                  min.variation,
                                  new = FALSE,
                                  param = 'PM2.5') {
  method <- match.arg(method)
  # Argument checks
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), length(data.column) == 1, data.column %in% colnames(data))
  stopifnot(persist.duration > 0, !is.null(persist.duration))
  stopifnot(compl.persist.duration >= 0, !is.null(compl.persist.duration))
  stopifnot(min.variation >= 0)
  stopifnot(method %in% c('range','std'))
  #
  compl.persist.duration <- ifelse(compl.persist.duration < 1,
                                   round(compl.persist.duration * persist.duration),
                                   compl.persist.duration)
  #
  constant_value_fun <- function(df) {
    ind.nona <- which(!is.na(df[[data.column]]))
    if (length(ind.nona) == 0) {
      df <- df %>%
        mutate(flag_constant_value = NA_real_)
    } else {
      start.nona <- df[[date.column]][ind.nona[1]]
      end.nona   <- df[[date.column]][ind.nona[length(ind.nona)]]
      df.sub <- df %>%
        filter(.data[[date.column]] >= start.nona, .data[[date.column]] <= end.nona)
      nobs.statistic <- zoo::rollapply(df.sub[[data.column]], persist.duration,
                                       FUN = function(x) sum(!is.na(x)),
                                       by.column = FALSE, fill = NA_real_, align = 'right')
      statistic <- if (method == "range") {
        rollapply(df.sub[[data.column]], persist.duration,
                       FUN = function(x) {
                         if (all(is.na(x))) return(NA_real_)
                         round(max(x, na.rm = TRUE) - min(x, na.rm = TRUE), 1)
                       },
                       by.column = FALSE, fill = NA_real_, align = 'right')
      } else {
        rollapply(df.sub[[data.column]], persist.duration,
                       FUN = function(x) round(sd(x, na.rm = TRUE), 1),
                       by.column = FALSE, fill = NA_real_, align = 'right')
      }
      df.sub <- df.sub %>%
        mutate(flag_constant_value = case_when(
          nobs.statistic < compl.persist.duration ~ 0,
          statistic <= min.variation & nobs.statistic >= compl.persist.duration ~ 1,
          statistic >  min.variation & nobs.statistic >= compl.persist.duration ~ 2,
          TRUE ~ NA_real_
        ))
      df.out <- df %>%
        filter(.data[[date.column]] < start.nona | .data[[date.column]] > end.nona)
      df <- bind_rows(df.out, df.sub) %>%
        arrange(.data[[date.column]]) %>%
        mutate(flag_constant_value = if_else(is.na(.data[[data.column]]), NA_real_, .data[["flag_constant_value"]]))
      }
    if (new) {
      df <- df %>%
        mutate(!!paste0(param, "_constant_value") := ifelse(.data[["flag_constant_value"]] %in% c(0, 1) | is.na(.data[["flag_constant_value"]]), NA_real_,.data[[data.column]]))
    }
    return(df)
  }
  # Run on groups or entire data
  if (!is.null(id.column)) {
    stopifnot(is.character(id.column), id.column %in% colnames(data))
    data <- data %>%
      group_by(across(all_of(id.column))) %>%
      group_map(~constant_value_fun(.x), .keep = TRUE) %>%
      bind_rows()
  } else {
    data <- constant_value_fun(data)
  }
  return(as_tibble(data))
}
