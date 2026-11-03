#' Total Quality Flag
#'
#' Calculates the overall quality flag for each measurement based on the individual quality flags from the different statistical tests.
#'
#' @details
#' The quality flag, \eqn{Q}, is determined as a function of the individual flags derived from the quality control algorithms.
#'
#' \deqn{
#' Q = h(f^\text{range}, f^\text{constant}, f^\text{outlier}, f^\text{correlation}, f^\text{similarity})
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{f^\text{range} \in \{ 0, 1 \}}
#'   \item \eqn{f^\text{constant} \in \{ 0, 1, 2 \}}
#'   \item \eqn{f^\text{outlier} \in \{ 0, 1, 2 \}}
#'   \item \eqn{f^\text{correlation} \in \{ 0, 1, 2, 3 \}}
#'   \item \eqn{f^\text{similarity} \in \{ 0, 1, 2, 3 \}}
#' }
#'
#' \eqn{Q}, is categorized into three levels: Other, Good, and High, based on the individual quality control flags as follows:
#'
#' \deqn{
#' Q = \begin{cases}
#' \text{High Quality} & (f^\text{range}, f^\text{constant}, f^\text{outlier}, f^\text{correlation}, f^\text{similarity}) = (1, 2, 2, 3, 3)\\
#' \text{Good Quality} & (f^\text{range}, f^\text{constant}, f^\text{outlier}, f^\text{correlation}) = (1, 2, 2, 3) \text{, } f^\text{similarity} \in \{ 0, 1, 2 \}\\
#' \text{Other Quality} & \text{otherwise}.
#' \end{cases}
#' }
#'
#'
#' @param data A data frame containing the dataset to be tested
#' @param data.column character. Name of the column in `data` on which the function is performed.
#' @param qcflags character. A vector including the variable names of the individual quality flags. Default is `c('flag_range','flag_constant_value', 'flag_outlier','flag_correlation','flag_similarity')`
#'
#' @return Returns the original dataset with:
#' \itemize{
#'   \item `quality_flag`: integer Total quality flag (`NA` - not available, `0` - other quality, `1` - good quality, and `2` - high quality)
#' }
#'
#' @seealso \code{\link{range_validity_qctest}}, \code{\link{constant_value_qctest}}, \code{\link{outlier_qctest}}, \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, and \code{\link{spatial_similarity_qctest}}.
#'
#'
#' @export
quality_flag_function <- function(data, data.column,
                                  qcflags = c("flag_range", "flag_constant_value","flag_outlier",
                                              "flag_correlation","flag_similarity")) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(data.column), length(data.column) == 1, data.column %in% colnames(data))
  stopifnot(is.character(qcflags) || is.null(qcflags))
  #
  use.similarity  <- ("flag_similarity"  %in% qcflags) && ("flag_similarity"  %in% colnames(data))
  use.correlation <- ("flag_correlation" %in% qcflags) && ("flag_correlation" %in% colnames(data))
  data <- data %>%
    mutate(quality_flag = case_when(
      is.na(.data[[data.column]]) ~ NA_integer_,
      !is.na(.data[[data.column]]) &
        (is.na(flag_range) | is.na(flag_constant_value) | is.na(flag_outlier)) ~ 0L,
      !is.na(.data[[data.column]]) & !("flag_correlation" %in% qcflags) ~ 0L,
      use.similarity &
        flag_range == 1 & flag_constant_value == 2 & flag_outlier == 2 &
        flag_correlation == 3 & (is.na(flag_similarity) | flag_similarity %in% c(0,1,2)) ~ 1L,
      use.similarity &
        flag_range == 1 & flag_constant_value == 2 & flag_outlier == 2 &
        flag_correlation == 3 & flag_similarity == 3 ~ 2L,
      !use.similarity &
        flag_range == 1 & flag_constant_value == 2 & flag_outlier == 2 &
        flag_correlation == 3 ~ 1L,
      TRUE ~ 0L
    ))
  return(as_tibble(data))
}
