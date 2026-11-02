#' Spatial Outlier Detection Test
#'
#' Verifies whether measurements flagged as outliers by \code{\link{outlier_qctest}} are actual outliers or represent a special event (e.g., a pollution episode), by checking whether neighboring locations show similar abnormal values at the same timestamps.
#' It updates the quality flag assigned by \code{\link{outlier_qctest}} and optionally creates a column of valid values.
#'
#'
#' @details
#' If a measurement is flagged as an outlier at a target location and neighboring locations also show outliers at the same timestamp, the measurement is not considered abnormal, since it may reflect an real event. In this case, the quality flag of \code{\link{outlier_qctest}} is updated for the target location.
#'
#' Let \eqn{X_{it}} and \eqn{f_{it}^{\text{outlier}}} denote the measurements and the quality flag from the \code{\link{outlier_qctest}} at a location \eqn{i} at timestamp \eqn{t}. The test is applied only at the timestamps flagged as outliers:
#' \deqn{
#' t^{\ast} = \{ t \mid f_{it}^{\text{outlier}} = 1 \}.
#' }
#'
#' At each timestamp \eqn{t^{\ast}}, the neighboring locations \eqn{j} within a radius of influence \eqn{\rho} are identified, and the number of neighbors that are also flagged as outliers is calculated as:
#' \deqn{
#' N_{\rho}(i,t^{\ast}) = \left| \{ j : d_{ij} \leq \rho \text{ and } f_{jt^{\ast}}^{\text{outlier}} = 1 \} \right|,
#' }
#' where \eqn{d_{ij}} is the separation distance between locations \eqn{i} and \eqn{j}.
#'
#' The quality flag \eqn{f_{it^{\ast}}^{\text{outlier}}} is then updated as:
#' \deqn{
#' f_{it^{\ast}}^{\text{outlier}} =
#' \begin{cases}
#' 1 & \text{Outlier: } N_{\rho_{\text{local}}}(i,t^{\ast}) < 1 \text{ and } N_{\rho_{\text{regional}}}(i,t^{\ast}) < n_{\text{min,neigh}}, \\
#' 2 & \text{Non-outlier: } N_{\rho_{\text{local}}}(i,t^{\ast}) \geq 1 \text{ or } N_{\rho_{\text{regional}}}(i,t^{\ast}) \geq n_{\text{min,neigh}},
#' \end{cases}
#' }
#' \itemize{
#'   \item \eqn{\rho_{\text{local}}}{rho_local}, \eqn{\rho_{\text{regional}}}{rho_regional} (`outlier.radii`): inner and outer radii of influence, in km.
#'   \item \eqn{n_{\text{min,neigh}}}{n_min,neigh} (`outlier.minneigh`): minimum number of neighbors at a distance \eqn{d_{ij}} with \eqn{\rho_{\text{local}} < d_{ij} \leq \rho_{\text{regional}}}{rho_local < d_ij <= rho_regional} required to consider the measurement a non-outlier.
#' }
#'
#' @param data A data frame containing the dataset to be tested
#' @param date.column character. Name of the date/time column in `data`.
#' @param data.column character. Name of the column in `data` to test.
#' @param out.flag.column character. Name of the outlier flag column in `data`, assigned by \code{\link{outlier_qctest}}.
#' @param id.column character. Name of the column identifying the measurement location ID.
#' @param x.column character. Name of the longitude column in `data`.
#' @param y.column character. Name of the latitude column in `data`.
#' @param outlier.radii numeric. Vector of two values with the inner (local) and the outer (regional) radius of influence, in km. Default is `c(2, 25)`.
#' @param outlier.minneigh numeric. Minimum number of neighbors outside the inner radius required to update the outlier flag. Default is `2`.
#' @param new logical. If TRUE, adds a column `<param>_outlier` containing the valid values; flagged values are replaced with `NA`. Default is `FALSE`.
#' @param param character. Base name of the new column (e.g. `PM2.5` gives `PM2.5_outlier`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_outlier}: integer. Quality flags based on the statistical test, \code{0} = insufficient data, \code{1} = outlier, \code{2} = non-outlier. It is `NA` where the tested value is missing.
#'   \item \code{<param>_outlier}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
spatial_outlier_qctest <- function(data,
                                   date.column,
                                   data.column,
                                   out.flag.column,
                                   id.column,
                                   x.column,
                                   y.column,
                                   outlier.radii = c(2, 25),
                                   outlier.minneigh = 2,
                                   new = FALSE,
                                   param = 'PM2.5') {
  # Argument checks
  stopifnot(all(sapply(c(date.column, data.column, out.flag.column, id.column, x.column, y.column), is.character)))
  stopifnot(all(c(date.column, data.column, out.flag.column, id.column, x.column, y.column) %in% colnames(data)))
  stopifnot(is.numeric(outlier.radii), length(outlier.radii) > 1, !is.null(outlier.radii))
  stopifnot(is.numeric(outlier.minneigh), !is.null(outlier.minneigh))
  stopifnot(!missing(param) && !is.null(param))
  # maximum outlier.radii in m
  outlier.radii <- outlier.radii * 1000
  #
  data <- data %>%
    group_by(across(all_of(id.column))) %>%
    group_map(~{
      id.target <- .y[[1]]
      df <- .x
      df.target <- df %>%
        filter(.data[[out.flag.column]] == 1)
      #
      if (nrow(df.target) > 0) {
        info.target <- df.target %>%
          select(all_of(c(id.column, x.column, y.column))) %>%
          distinct()
        sp.target <- st_as_sf(info.target, coords = c(x.column, y.column), crs = 4326L)
        # neighboring locations
        # neighbor data
        timestamps <- unique(df.target[[date.column]])
        df.neigh <- data %>%
          filter(.data[[id.column]] != id.target,
                 .data[[date.column]] %in% timestamps,
                 .data[[out.flag.column]] == 1)
        #
        if (nrow(df.neigh) > 0) {
          info.neigh <- data %>%
            filter(.data[[id.column]] %in% unique(df.neigh[[id.column]])) %>%
            select(all_of(c(id.column, x.column, y.column))) %>%
            distinct()
          sp.neigh <- st_as_sf(info.neigh, coords = c(x.column, y.column), crs = 4326L)
          # calculation of separation distance and selection of neighbors
          distances <- as.vector(round(st_distance(sp.target, sp.neigh), 1))
          info.neigh <- info.neigh %>%
            mutate(dij = distances) %>%
            filter(dij <= outlier.radii[2])
          #
          if (nrow(info.neigh) > 0) {
            # calculation of the spatial outlier condition
            # local scale: <= 2 km
            # not-local scale: > 2 km
            df.neigh <- df.neigh %>%
              filter(.data[[id.column]] %in% info.neigh[[id.column]])
            df.neigh <- left_join(df.neigh, info.neigh %>% select(all_of(c(id.column, "dij"))), by = id.column)
            df.target.neigh <- left_join(df.target %>% select(all_of(c(date.column, id.column))) %>% rename(!!paste0(id.column, "_target") := all_of(id.column)),
                                         df.neigh %>% select(all_of(c(date.column, id.column, "dij"))) %>% rename(!!paste0(id.column, "_neigh") := all_of(id.column)),
                                         by = date.column)
            #
            if (nrow(df.target.neigh) > 0) {
              df.target.neigh <- df.target.neigh %>%
                mutate(radii_cut = ifelse(dij <= outlier.radii[1], 'local', 'not_local'))
              spatial.stat <- df.target.neigh %>%
                group_by(!!sym(date.column)) %>%
                summarise(local = length(which(radii_cut == 'local')),
                          not_local = length(which(radii_cut == 'not_local'))) %>%
                mutate(!!paste0(out.flag.column, "_spatial") := ifelse(local >= 1 | not_local >= outlier.minneigh, 2, 1))
              #
              df <- df %>%
                mutate(!!out.flag.column := if_else(.data[[date.column]] %in% spatial.stat[[date.column]],
                                                    spatial.stat[[paste0(out.flag.column, "_spatial")]][match(.data[[date.column]], spatial.stat[[date.column]])], .data[[out.flag.column]]))
            }
          }
        }
      }
      # Create {param}_outlier column once at the end
      if (new == TRUE) {
        df <- df %>%
          mutate(!!paste0(param, '_outlier') := ifelse(.data[[out.flag.column]] == 0 | .data[[out.flag.column]] == 1 | is.na(.data[[out.flag.column]]), NA_real_, .data[[data.column]]))
      }
      return(as_tibble(df))
    }, .keep = TRUE) %>%
    bind_rows()
}
