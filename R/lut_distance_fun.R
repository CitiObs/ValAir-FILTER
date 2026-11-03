#' Correction Distance by Season
#'
#' Internal helper function that returns the maximum distance (in km) up to which the correction is performed, according to the season of the year.
#'
#' @details
#' The distance is retrieved from a look-up table and corresponds to the season-dependent radius of influence \eqn{\rho_{\text{corr}}}{rho_corr} used by \code{\link{correction_process}} to select the reference stations. Seasons are coded as follows:
#'
#' 1. Winter: December 1 to February 28 (29 in a leap year)
#' 2. Spring: March 1 to May 31
#' 3. Summer: June 1 to August 31
#' 4. Autumn (Fall): September 1 to November 30
#'
#' @param season integer. Season of the year, coded from `1` to `4` as described in Details.
#'
#' @return A numeric value: the distance (in km) for the specified season.
#'
#' @export
lut_distance_fun <- function(season) {
  Season = c(1, 2, 3, 4)
  Distance = c(11.4847, 12.71141, 19.9706, 17.08543)  # in km
  lookup_distance <- tibble(Season = Season, Distance = Distance)
  dist <- (lookup_distance %>% dplyr::filter(Season %in% season))[['Distance']]
  return(dist)
}
