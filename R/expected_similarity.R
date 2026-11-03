#' Expected Similarity
#'
#' Estimates the expected Euclidean distance between PM\eqn{_{2.5}} concentrations at two locations over a week as a function of the geographic distance between the two locations and calendar season.
#' It is calculated using the PM\eqn{_{2.5}} concentrations from the European Environment Agency (EEA) air quality monitoring sites distributed across Europe.
#'
#' @details
#' For a given season \eqn{s} and distance \eqn{d}, the expected Euclidean distance \eqn{D_{exp}(d,s)} is computed as:
#'
#' \deqn{D_{exp}(d,s) = b_0(s) + b_1(s) d + b_2(s) d^2 + b_3(s) d^3}
#'
#' where \eqn{s \in \{1, 2, 3, 4\}} represents the seasons
#'
#' 1. Winter: December 1 to February 28 (29 in a leap year)
#' 2. Spring: March 1 to May 31
#' 3. Summer: June 1 to August 31
#' 4. Autumn (Fall): September 1 to November 30
##'
#' and the season-specific coefficients \eqn{b_0(s), b_1(s), b_2(s), b_3(s)} are given as:
#'
#' \deqn{
#' b_0 = (68, 58.5, 50.7, 57.1)
#' }
#' \deqn{
#' b_1 = (0.502, 0.453, 0.286, 0.333)
#' }
#' \deqn{
#' b_2 = (-1.3 \times 10^{-3}, -1.26 \times 10^{-3}, -0.85 \times 10^{-3}, -0.91 \times 10^{-3})
#' }
#' \deqn{
#' b_3 = (1.242 \times 10^{-6}, 1.213 \times 10^{-6}, 8.415 \times 10^{-7}, 8.788 \times 10^{-7})
#' }
#'
#' A relaxation factor \eqn{r} can also be applied, and the expected Euclidean distance \eqn{D_{exp}(d,s)} is adjusted as:
#'
#' \deqn{
#' D_{exp}(d,s) = D_{exp}(d,s) \times (1 - r)
#' }
#'
#' @param season integer. Season of the year (`1` = Winter, `2` = Spring, `3` = Summer, and `4` = Autumn).
#' @param distance numeric. A numeric value specifying the separation distance (in meters, m).
#' @param relax.factor numeric. Default is `0.25`.
#'
#' @return Returns expected similarity values
#'
#'
#' @export
expected_similarity <- function(season, distance, relax.factor = 0.25) {
  Season <- c(1,2,3,4)
  b0 <- c(68.04585477, 58.51618754, 50.71076079, 57.14827933)
  b1 <- c(0.50238293, 0.452549077, 0.286315948, 0.333297642)
  b2 <- c(-0.001295603, -0.001263275, -0.000851074, -0.000914356)
  b3 <- c(1.24202946478091e-06, 1.21303177254065e-06, 8.4149769556485e-07, 8.78846013819774e-07)
  lookup_similarity_distance <- tibble(Season = Season, b0, b1, b2, b3)
  #
  params <- tibble(Season = season, distance = distance)
  params <- left_join(params, lookup_similarity_distance, by = 'Season') %>% as_tibble()
  similarity <- params[['b0']] + params[['b1']]*params[['distance']] + params[['b2']]*params[['distance']]^2 + params[['b3']]*params[['distance']]^3
  if (!is.null(relax.factor)) {
    similarity <- similarity - relax.factor*similarity
  }
  return(similarity)
}
