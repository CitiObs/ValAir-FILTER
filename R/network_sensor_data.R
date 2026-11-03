#' PM\eqn{_{2.5}} Observations from a Low-cost Sensor Network
#'
#' Sample hourly PM\eqn{_{2.5}} concentrations from a network of low-cost sensors.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{\code{id}}{numeric. ID of each measurement location in the network.}
#'   \item{\code{timestamp}}{POSIXct. Hourly timestamps in \code{"YYYY-MM-DD HH:MM:SS"} format. Timestamps are in UTC.}
#'   \item{\code{longitude}}{numeric. Longitude of each measurement location (decimal degrees).}
#'   \item{\code{latitude}}{numeric. Latitude of each measurement location (decimal degrees).}
#'   \item{\code{PM2.5}}{numeric. Hourly PM\eqn{_{2.5}} concentration values (\eqn{\mu g \, m^{-3}}).}
#' }
#'
#' @usage data(network_sensor_data)
#' @keywords datasets
#'
#' @source Sample data for demonstration purposes.
'network_sensor_data'