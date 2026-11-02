#' Inverse Distance Weighted Interpolation
#'
#' Calculates the weighted values of one or multiple parameters from neighboring locations, using the Inverse Distance Weighted (IDW) interpolation with spatial weights based on the modified Shepard method.
#'
#'
#' @details
#' Let \eqn{X_{jt}} denote the measurements at each timestamp \eqn{t \in \mathcal{T}} of the locations \eqn{j \in \mathcal{J}}, which are the neighbors of the measurement location \eqn{i}. Missing values are excluded from the calculation. The weighted time series \eqn{X_{wt}} is calculated as:
#' \deqn{
#' X_{wt} = \frac{\sum_{j \in \mathcal{J}} w_{ij} X_{jt}}{\sum_{j \in \mathcal{J}} w_{ij}},
#' }{X_w = sum(w_ij * X_jt) / sum(w_ij)}
#' where the spatial weight \eqn{w_{ij}} is defined by the modified Shepard method as:
#' \deqn{
#' w_{ij} =
#' \begin{cases}
#' \left( \frac{\rho - d_{ij}}{\rho \, d_{ij}} \right)^2, & d_{ij} < \rho, \\
#' 0, & d_{ij} \geq \rho,
#' \end{cases}
#' }{w_ij = ((rho - d_ij) / (rho * d_ij))^2 if d_ij < rho, 0 otherwise}
#' \itemize{
#'   \item \eqn{w_{ij}}: spatial weight assigned to the locations \eqn{i} and \eqn{j}.
#'   \item \eqn{d_{ij}}: separation distance between \eqn{i} and \eqn{j}.
#'   \item \eqn{\rho} (`radii`): radius of influence, in km.
#' }
#'
#' @param info A data frame containing the spatial information (e.g., ID, longitude, latitude, separation distance) of the measurement location, e.g., the output of \code{\link{find_neighbors}}.
#' @param data A data frame containing the dataset of the measurement locations, e.g., the output of \code{\link{find_neighbors}}.
#' @param radii numeric. Radius of influence, in km, used to calculate the weights in the modified Shepard weighting function.
#' @param date.column character. Name of the date/time column in `data`.
#' @param id.column character. Name of the column identifying the measurement location ID in `data`.
#' @param data.column character. Name of the column(s) in `data` on which the function is performed.
#'
#' @return A data frame with the date/time column and the interpolated values.
#'
#' @seealso \code{\link{find_neighbors}} and \code{\link{correction_process}}.
#'
#' @export
modified_idw_function <- function(info,
                                  data,
                                  radii,
                                  date.column,
                                  id.column,
                                  data.column) {
  # argument flags
  stopifnot(!is.null(info),!is.null(data))
  stopifnot(is.numeric(radii), !is.null(radii))
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(id.column), id.column %in% colnames(data))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  # Distance matrix and weights
  distance <- matrix(info[['dij']], nrow = 1)
  w <- ((1000 * radii - distance) / (1000 * radii * distance))^2
  for (i in 1:nrow(w)) {
    if (sum(is.infinite(w[i, ])) > 0) {
      w[i, !is.infinite(w[i, ])] <- 0
      w[i, is.infinite(w[i, ])] <- 1
    }
  }
  data.out <- data %>%
    dplyr::select(all_of(date.column)) %>%
    distinct()
  for (i in 1:length(data.column)) {
    mat.data <- pivot_wider(
      data %>% select(all_of(c(date.column, id.column, data.column[i]))),
      names_from = all_of(id.column),
      values_from = all_of(data.column[i])) %>% as_tibble() %>%
      select(-all_of(date.column))
    mat.w <- matrix(rep(as.vector(w),each=nrow(mat.data)), nrow=nrow(mat.data))
    ind.mat.w <- apply(mat.data, 2, function(x) ifelse(is.na(x), FALSE, TRUE))
    mat.w[!ind.mat.w] <- NA
    variable.weighted <- round(apply(mat.w * mat.data, 1, sum, na.rm = TRUE)/apply(mat.w, 1, sum, na.rm = TRUE),1)
    data.out <- data.out %>%
      mutate(!!data.column[i] := variable.weighted)
  }
  return(as_tibble(data.out))
}
