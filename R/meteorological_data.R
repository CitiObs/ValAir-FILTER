#' Meteorological Observations from Weather Stations
#'
#' Sample hourly air temperature and relative humidity from multiple NOAA meteorological stations.
#' 
#' @format A data frame with seven columns:
#' \describe{
#'   \item{\code{code}}{character. The ID of each station.}
#'   \item{\code{timestamp}}{POSIXct. Hourly timestamps in \code{"YYYY-MM-DD HH:MM:SS"} format. Timestamps are in UTC.}
#'   \item{\code{longitude}}{numeric. Longitude of each station (decimal degrees).}
#'   \item{\code{latitude}}{numeric. Latitude of each station (decimal degrees).}
#'   \item{\code{Temperature}}{numeric. Hourly air temperature observations (\eqn{^\circ}C)}
#'   \item{\code{Humidity}}{numeric. Hourly relative humidity observations (%).}
#' }
#'
#' @usage data(meteorological_data)
#' @keywords datasets
#'
#' @source Sample data for demonstration purposes
'meteorological_data'
