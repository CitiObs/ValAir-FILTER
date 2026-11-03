#' PM\eqn{_{2.5}} Observations from Reference Stations
#'
#' Sample hourly PM\eqn{_{2.5}} concentrations from multiple official monitoring stations.
#' 
#' @format A data frame with six columns:
#' \describe{
#'   \item{\code{id}}{integer. The ID of each station.}
#'   \item{\code{timestamp}}{POSIXct. Hourly timestamps in \code{"YYYY-MM-DD HH:MM:SS"} format. Timestamps are in UTC.}
#'   \item{\code{longitude}}{numeric. Longitude of each station (decimal degrees).}
#'   \item{\code{latitude}}{numeric. Latitude of each station (decimal degrees).}
#'   \item{\code{value}}{numeric. Hourly PM\eqn{_{2.5}} concentration values (\eqn{\mu g \, m^{-3}}).}
#' }
#'
#' @usage data(reference_data)
#' @keywords datasets
#'
#' @source Simulated data for demonstration purposes
'reference_data'
