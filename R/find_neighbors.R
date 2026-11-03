#' Find Nearest Neighbors
#'
#' Extracts the nearest neighbors of a specific location up to a specified distance, provided that they have enough data.
#'
#'
#' @details
#' Let \eqn{\mathcal{J}} be the set of the neighboring measurement locations \eqn{j} of a site \eqn{i}. The selection is performed in two steps:
#'
#' \strong{Step 1: Validity condition}
#'
#' Neighbors are required to satisfy both the data completeness and the distance constraints:
#' \deqn{
#' \mathcal{J}^{\ast} = \left\{ j \in \mathcal{J} \mid N_{\mathcal{T}_j} \geq n_{\text{min}} \text{ and } d_{ij} \leq \rho \right\},
#' }{J* = j in J with N_j >= n_min and d_ij <= rho}
#' \itemize{
#'   \item \eqn{\mathcal{J}^{\ast}}: subset of \eqn{\mathcal{J}} including only the locations \eqn{j} that satisfy both criteria.
#'   \item \eqn{N_{\mathcal{T}_j}}: number of available (non-missing) measurements of \eqn{j} within the time range \eqn{\mathcal{T}}.
#'   \item \eqn{n_{\text{min}}}{n_min} (`compl.threshold`): minimum number of non-missing measurements required, given either as a number of measurements or as a fraction of the time range (values `< 1`).
#'   \item \eqn{d_{ij}}: geographic distance between \eqn{i} and \eqn{j}.
#'   \item \eqn{\rho} (`maxradii`): radius of influence, i.e., the maximum distance considered, in km.
#' }
#'
#' \strong{Step 2: Selection condition}
#'
#' If a maximum number of neighbors \eqn{n_{\text{max,neigh}}}{n_max,neigh} (`maxneigh`) is specified, only the closest neighbors are retained:
#' \deqn{
#' \mathcal{J}^{\ast} =
#' \begin{cases}
#' \mathcal{J}^{\ast}, & \text{if } n_{\text{max,neigh}} \text{ is not specified or } |\mathcal{J}^{\ast}| \leq n_{\text{max,neigh}}, \\
#' \{ j_{(1)}, j_{(2)}, \dots, j_{(n_{\text{max,neigh}})} \}, & \text{otherwise},
#' \end{cases}
#' }{J* = the n_max,neigh closest locations of J*}
#' where \eqn{j_{(k)}} is the \eqn{k}-th closest neighbor of site \eqn{i}, ordered by increasing distance \eqn{d_{ij}}.
#'
#' @param data A data frame containing the dataset of the target locations.
#' @param date.column character. Name of the date/time column in `data` and `rdata`.
#' @param data.column character. Name of the column in `data` on which the function is performed.
#' @param id.column character. Name of the column identifying the measurement location ID in `data`.
#' @param rdata A data frame containing the dataset of the additional locations (e.g., reference stations).
#' @param rdata.column character. Name of the column in `rdata` used in the function.
#' @param rid.column character. Name of the column identifying the measurement location ID in `rdata`.
#' @param x.column character. Name of the longitude column in `data` and `rdata`.
#' @param y.column character. Name of the latitude column in `data` and `rdata`.
#' @param compl.threshold numeric. Minimum data completeness required to select a neighbor, expressed either as a number of measurements or as a fraction of the time window. Default is `90`.
#' @param maxradii numeric. Radius of influence, in km, used to extract the nearest neighbors. Default is `30`.
#' @param maxneigh numeric. Maximum number of neighboring locations. Default is `NULL`.
#'
#' @return A list of data frames:
#' \itemize{
#'   \item \code{data}: subset of `rdata` with the measurements of the selected neighboring locations.
#'   \item \code{neighbors}: data frame with the spatial information (e.g., ID, longitude, latitude) and the separation distance (`dij`) of the selected nearest neighbors.
#' }
#'
#' @seealso \code{\link{correction_process}} for the main calibration routine that uses this helper function.
#'
#' @export
find_neighbors <- function(data,
                           date.column,
                           data.column,
                           id.column,
                           rdata,
                           rdata.column,
                           rid.column,
                           x.column,
                           y.column,
                           maxradii,
                           maxneigh,
                           compl.threshold) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(data), date.column %in% colnames(rdata))
  stopifnot(is.character(id.column), id.column %in% colnames(data))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(x.column),x.column %in% colnames(data), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column),y.column %in% colnames(data),  y.column %in% colnames(rdata))
  stopifnot(is.numeric(maxradii), !is.null(maxradii))
  stopifnot(is.numeric(compl.threshold), !is.null(compl.threshold))
  #
  if (id.column == rid.column) {
    rdata <- rdata %>%
      rename(!!paste0('r_',rid.column) := rid.column)
    rid.column <- paste0('r_',rid.column)
  }
  # target location
  info.target <- data %>%
    select(all_of(c(id.column, x.column, y.column))) %>% distinct()
  sp.target <- st_as_sf(info.target, crs = 4326L, coords = c(x.column,y.column))
  # neighbors reference data
  df.neigh <- rdata %>%
    filter(.data[[date.column]] %in% data[[date.column]])
  info.neigh <- df.neigh %>%
    select(all_of(c(rid.column, x.column, y.column))) %>% distinct()
  if (nrow(info.neigh) == 0) {
    return(list(data=NULL,neighbors=NULL))
  } else {
    sp.neigh <- st_as_sf(info.neigh, crs = 4326L, coords = c(x.column,y.column))
    # calculation of separation distance and selection of neighbors
    distances <- round(st_distance(sp.target, sp.neigh),1)
    info.neigh <- info.neigh %>%
      mutate(dij = as.vector(distances)/1000) %>%
      arrange(.data[['dij']]) %>%
      filter(.data[['dij']] > 0 & .data[['dij']] < maxradii)
    if (nrow(info.neigh) == 0) {
      return(list(data=NULL,neighbors=NULL))
    } else {
      info.neigh <- info.neigh %>%
        arrange(.data[['dij']])
      df.neigh <- df.neigh %>%
        filter(.data[[rid.column]] %in% info.neigh[[rid.column]])
      df.target.neigh <- left_join(df.neigh, data %>% select(all_of(c(date.column, data.column, id.column))), by=date.column)
      #
      # data completeness and selection of neighboring location
      df.completeness <- df.target.neigh %>%
        select(all_of(c(rid.column, rdata.column, data.column))) %>%
        drop_na(all_of(c(rdata.column, data.column))) %>%
        group_by(.data[[rid.column]]) %>%
        summarise(across(everything(), ~ sum(!is.na(.x), na.rm = TRUE)), .groups = "drop")
      df.completeness <- df.completeness %>%
        filter(.data[[data.column]] >= compl.threshold)
      if (nrow(df.completeness) == 0) {
        return(list(data=NULL,neighbors=NULL))
      } else {
        df.completeness <- left_join(df.completeness, info.neigh %>% select(all_of(c(rid.column,'dij'))), by = rid.column) %>% as_tibble()
        df.completeness <- df.completeness %>%
          dplyr::arrange(.data[['dij']])
        if (is.null(maxneigh)) {
          maxneigh <- nrow(df.completeness)
        }
        # selection of the nearest neighboring locations
        closest.neigh <- df.completeness %>%
          arrange(.data[['dij']]) %>%
          slice_head(n = maxneigh)
        df.target.neigh <- df.target.neigh %>%
          filter(.data[[rid.column]] %in% closest.neigh[[rid.column]])
        info.target.neigh <- info.neigh %>%
          filter(.data[[rid.column]] %in% unique(df.target.neigh[[rid.column]]))
        #
        return(list(data=df.target.neigh,neighbors=info.target.neigh))
      }
    }
  }
}
