#' Range Validity Check
#'
#' Checks if measurements are within physically plausible limits, assigns a quality flag, and optionally creates a column of valid values.
#'
#' @details
#' Let \eqn{X_t} the measurement at timestamp \eqn{t}. The quality flag \eqn{f_{t}^{\mathrm{range}}} is defined as:
#' \deqn{
#'   f_{t}^{\mathrm{range}} =
#'   \begin{cases}
#'   0 & \text{Out of Range: } X_t < th_l \text{ or } X_t > th_u, \\
#'   1 & \text{Within Range: } th_l \leq X_t \leq th_u
#'   \end{cases}
#' }
#'
#' with:
#' \itemize{
#'   \item \eqn{th_l}: lower threshold (`lower.bound`)
#'   \item \eqn{th_u}: upper threshold (`upper.bound`)
#' }
#'
#'
#' @param data A data frame containing the dataset to be tested
#' @param data.column character. Name of the column in `data` to test.
#' @param upper.bound numeric. The upper threshold limit.
#' @param lower.bound numeric. The lower threshold limit.
#' @param new logical. If `TRUE`, creates a new column named `<param>_range` (e.g., `PM2.5_range`) that contains only the valid values; values outside the defined limits are replaced with `NA`. Default is `FALSE`.
#' @param param character. The base name to use for the `new` column (e.g., `PM2.5` becomes `PM2.5_range`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#'
#' @return Returns the original dataset with:
#' \itemize{
#'   \item \code{flag_range}: integer. Quality flags based on the statistical test, \code{0} = out-of-range, \code{1} = within range. It is `NA` where the tested value is missing.
#'   \item \code{<param>_range}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
range_validity_qctest <- function(data,
                                  data.column,
                                  upper.bound,
                                  lower.bound,
                                  new = FALSE,
                                  param = 'PM2.5') {
  # argument checks
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  stopifnot(is.numeric(upper.bound), is.numeric(lower.bound),
            !is.null(upper.bound), !is.null(lower.bound),
            upper.bound > lower.bound)
  # range test application
  data <- data %>%
    mutate(flag_range = ifelse(is.na(.data[[data.column]]), NA_real_,
                               ifelse(.data[[data.column]] >= lower.bound & .data[[data.column]] <= upper.bound,  1, 0)))
  if (new) {
    data <- data %>% mutate(!!paste0(param, "_range") := ifelse(.data[['flag_range']] == 0 | is.na(.data[['flag_range']]), NA_real_, .data[[data.column]]))
  }
  return(as_tibble(data))
}
