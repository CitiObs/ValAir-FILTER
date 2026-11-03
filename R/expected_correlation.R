#' Expected Correlation
#'
#' Estimates the expected Pearson's \eqn{R} correlation coefficient between PM\eqn{_{2.5}} concentrations at two locations over one month as a function of the geographic distance between the two locations and calendar season.
#' It is calculated using the PM\eqn{_{2.5}} concentrations from the European Environment Agency (EEA) air quality monitoring sites distributed across Europe.
#'
#' @details
#' For a given season \eqn{s} and distance \eqn{d}, the expected correlation coefficient \eqn{R_{exp}(d,s)} is computed as:
#'
#' \deqn{
#' R_{exp}(d,s) = b_0(s) + b_1(s) d + b_2(s) d^2 + b_3(s) d^3
#' }
#'
#' where \eqn{s \in \{1, 2, 3, 4\}} represents the seasons,
#'
#' 1. Winter: December 1 to February 28 (29 in a leap year)
#' 2. Spring: March 1 to May 31
#' 3. Summer: June 1 to August 31
#' 4. Autumn (Fall): September 1 to November 30
#'
#' and the season-specific coefficients \eqn{b_0(s), b_1(s), b_2(s), b_3(s)} are given as:
#'
#' \deqn{
#' b_0 = (0.867, 0.852, 0.764, 0.821)
#' }
#' \deqn{
#' b_1 = (-0.00265, -0.00301, -0.00326, -0.00286)
#' }
#' \deqn{
#' b_2 = (5.74 \times 10^{-6}, 6.26 \times 10^{-6}, 7.89 \times 10^{-6}, 6.34 \times 10^{-6})
#' }
#' \deqn{
#' b_3 = (-5.22 \times 10^{-9}, -5.11 \times 10^{-9}, -7.20 \times 10^{-9}, -5.53 \times 10^{-9})
#' }
#'
#' A relaxation factor \eqn{r} can also be applied, and the expected correlation, \eqn{R_{exp}(d,s)}, is adjusted as:
#'
#' \deqn{
#' R_{exp}(d,s) = R_{exp}(d,s) \times (1 - r)
#' }
#'
#'
#' @param season integer. Season of the year (`1` = Winter, `2` = Spring, `3` = Summer, and `4` = Autumn).
#' @param distance numeric. A numeric value specifying the separation distance (in meters, m).
#' @param relax.factor numeric. Default is `0.25`.
#'
#' @return Returns expected correlation values
#'
#' @export
expected_correlation <- function(season, distance, relax.factor = 0.25) {
  Season <- c(1,2,3,4)
  b0 <- c(0.866920787, 0.851493359, 0.764202409, 0.82053746)
  b1 <- c(-0.00265353,-0.00301283,-0.003258162,-0.00285821)
  b2 <- c(5.735157e-06, 6.257868e-06, 7.885344e-06, 6.338316e-06)
  b3 <- c(-5.216934e-09, -5.111353e-09, -7.204441e-09, -5.53291e-09)
  lookup_correlation_distance <- tibble(Season = Season, b0, b1, b2, b3)
  #
  params <- tibble(Season = season, distance = distance)
  params <- left_join(params, lookup_correlation_distance, by = 'Season') %>% as_tibble()
  correlation <- params[['b0']] + params[['b1']]*params[['distance']] + params[['b2']]*params[['distance']]^2 + params[['b3']]*params[['distance']]^3
  if (!is.null(relax.factor)) {
    correlation <- correlation - relax.factor*correlation
  }
  return(correlation)
}
