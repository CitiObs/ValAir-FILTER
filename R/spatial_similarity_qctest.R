#' Spatial Similarity Test
#'
#' Checks whether the measurements at a location are similar to those of neighboring official monitoring stations, in terms of the Euclidean distance, and assigns a quality flag. Optionally creates a column of valid values.
#'
#'
#' @details
#' Let \eqn{X_{it}} denote the measurements at a location \eqn{i}. For each time step
#' \eqn{t_i}, the temporal window \eqn{\mathcal{T}} is defined as:
#' \deqn{
#' \mathcal{T} = [t_i - D, t_i],
#' }{T = [t_i - D, t_i]}
#' where \eqn{D} (`similarity.duration`) is the duration of the window.
#'
#' Missing values are excluded from the calculation. Let \eqn{\mathcal{T}^{\ast} \subseteq \mathcal{T}} be the set of timestamps in the window with available
#' (non-missing) values for both \eqn{X_{it}} and the reference station measurements \eqn{Y_{jt}}, and \eqn{N_{\mathcal{T}^{\ast}}} their number. The Euclidean distance \eqn{D_{\mathcal{T}}(i,j)} between location \eqn{i} and each neighboring reference station \eqn{j}, together with their separation distance \eqn{d_{ij}}, is calculated as:
#' \deqn{
#' D_{\mathcal{T}}(i,j) = \sqrt{\sum\limits_{t \in \mathcal{T}^{\ast}} (X_{it} - Y_{jt})^2}.
#' }{D(i,j) = sqrt(sum((X_it - Y_jt)^2))}
#'
#' The calculated \eqn{D_{\mathcal{T}}(i,j)} is compared with an expected similarity \eqn{D_{\text{exp}}(d_{ij}, s)}, derived from a look-up table based on the separation distance \eqn{d_{ij}} and the season \eqn{s} of \eqn{t \in \mathcal{T}}. See \code{\link{expected_similarity}} for more information.
#'
#' For each window \eqn{\mathcal{T}}, the number of reference stations \eqn{j} within the radius of influence \eqn{\rho} (`similarity.radii`) that satisfy the following conditions is derived:
#' \itemize{
#'   \item Adequate pairwise measurements:
#'   \deqn{
#'   N_{j}(\rho) = \left| \{ j : N_{\mathcal{T}^{\ast}} \geq n_{\text{min}} \text{ and } d_{ij} \leq \rho \} \right|
#'   }{N_j = number of j with N_T* >= n_min and d_ij <= rho}
#'   \item Similarity below the expected threshold, \eqn{D_{\mathcal{T}}(i,j) < D_{\text{exp}}(d_{ij}, s)}:
#'   \deqn{
#'   N_{D_{\mathcal{T}}(i,j) < D_{\text{exp}}}(\rho) = \left| \{ j : D_{\mathcal{T}}(i,j) < D_{\text{exp}}(d_{ij}, s) \text{ and } d_{ij} \leq \rho \} \right|
#'   }{N_sim = number of j with D(i,j) < D_exp and d_ij <= rho}
#' }
#'
#' The quality flag \eqn{f_{it}^{\text{similarity}}} for the location \eqn{i} at each timestamp \eqn{t} is then defined as:
#' \deqn{
#' f_{it}^{\text{similarity}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: } N_{\mathcal{T}^{\ast}} < n_{\text{min}}, \\
#' 1 & \text{Isolated site: } N_{j}(\rho) < 1, \\
#' 2 & \text{No similarity: } N_{D_{\mathcal{T}}(i,j) < D_{\text{exp}}}(\rho) < 1, \\
#' 3 & \text{Similarity: } N_{D_{\mathcal{T}}(i,j) < D_{\text{exp}}}(\rho) \geq 1,
#' \end{cases}
#' }{f = 0 (insufficient data), 1 (isolated site), 2 (no similarity), 3 (similarity)}
#' \itemize{
#'   \item \eqn{n_{\text{min}}}{n_min} (`compl.similarity.threshold`): minimum number of non-missing values required to evaluate the test, given either as a number of measurements or as a fraction \eqn{p} of the window duration \eqn{D} (`similarity.duration`), in which case \eqn{n_{\text{min}} = p D}{n_min = p * D}.
#'   \item \eqn{\rho} (`similarity.radii`): radius of influence, in km.
#' }
#'
#' @param data A data frame containing the dataset to be tested
#' @param date.column character. Name of the date/time column in `data` and `rdata`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location ID in `data`. If `NULL`, `data` is treated as a single location. Default is `NULL`.
#' @param rdata A data frame containing the data of the reference stations.
#' @param rdata.column character. Name of the column in `rdata` used in the test.
#' @param rid.column character. Name of the column identifying the measurement location in `rdata`.
#' @param x.column character. Name of the longitude column in `data` and `rdata`.
#' @param y.column character. Name of the latitude column in `data` and `rdata`.
#' @param similarity.radii numeric. Radius of influence, in km. Default is `30`.
#' @param similarity.duration numeric. Length of the rolling window. Default is `168`.
#' @param compl.similarity.threshold numeric. Minimum data completeness required within `similarity.duration` to run the test, expressed either as a number of measurements or as a fraction of `similarity.duration`. Default is `90`.
#' @param new logical. If `TRUE`, adds a column `<param>_similarity` containing the valid values; flagged values are replaced with `NA`. Default is `FALSE`.
#' @param param character. Base name of the new column (e.g. `PM2.5` gives `PM2.5_similarity`). Default is `PM2.5`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_similarity}: integer. Quality flags based on the statistical test, \code{0} = insufficient data, \code{1} = isolated site, \code{2} = no similarity with the neighboring reference stations, \code{3} = similarity with the neighboring reference stations. It is `NA` where the tested value is missing.
#'   \item \code{<param>_similarity}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
spatial_similarity_qctest <- function(data,
                                      date.column,
                                      data.column,
                                      id.column = NULL,
                                      rdata,
                                      rdata.column,
                                      rid.column,
                                      x.column,
                                      y.column,
                                      similarity.radii = 30,
                                      similarity.duration = 168,
                                      compl.similarity.threshold = 90,
                                      new = FALSE,
                                      param = 'PM2.5') {
  # Early exit if rdata is NULL
  if (is.null(rdata)) {
    data <- data %>%
      mutate(flag_similarity = NA_real_,
             !!paste0(param, "_similarity") := NA_real_)
    return(data)
  }
  # Argument checks
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(x.column), x.column %in% colnames(data), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column), y.column %in% colnames(data), y.column %in% colnames(rdata))
  stopifnot(is.numeric(similarity.radii))
  stopifnot(is.numeric(similarity.duration), similarity.duration > 0)
  stopifnot(is.numeric(compl.similarity.threshold), compl.similarity.threshold > 0)
  stopifnot(!missing(param) && !is.null(param))
  if (!is.null(id.column)) {
    stopifnot(is.character(id.column), id.column %in% colnames(data))
  }
  #
  compl.similarity.threshold <- ifelse(compl.similarity.threshold < 1,
                                       round(compl.similarity.threshold * similarity.duration),
                                       compl.similarity.threshold)
  similarity.radii <- similarity.radii * 1000
  #
  data <- data %>%
    mutate(!!id.column := as.character(.data[[id.column]]))
  rdata <- rdata %>%
    mutate(!!rid.column := as.character(.data[[rid.column]]))
  #
  similarity_fun <- function(df) {
    #
    info.target <- df %>%
      select(all_of(c(if (!is.null(id.column)) id.column else NULL, x.column, y.column))) %>%
      distinct()
    sp.target <- sf::st_as_sf(info.target, coords = c(x.column, y.column), crs = 4326)
    ind.nona <- which(!is.na(df[[data.column]]))
    if (length(ind.nona) == 0) {
      df <- df %>% mutate(flag_similarity = NA_real_)
    } else {
      start.nona <- df[[date.column]][ind.nona[1]]
      end.nona <- df[[date.column]][ind.nona[length(ind.nona)]]
      df.sub <- df %>%
        filter(.data[[date.column]] >= start.nona, .data[[date.column]] <= end.nona)
      if (sum(!is.na(df.sub[[data.column]])) < compl.similarity.threshold) {
        df <- df %>% mutate(flag_similarity = 0)
      } else {
        df.neigh <- rdata %>%
          filter(.data[[date.column]] >= start.nona, .data[[date.column]] <= end.nona)
        if (nrow(df.neigh) == 0) {
          df <- df %>% mutate(flag_similarity = 1)
        } else {
          info.neigh <- df.neigh %>%
            select(all_of(c(rid.column, x.column, y.column))) %>%
            distinct()
          sp.neigh <- st_as_sf(info.neigh, coords = c(x.column, y.column), crs = 4326)
          distances <- as.numeric(st_distance(sp.target, sp.neigh))
          info.neigh <- info.neigh %>%
            mutate(dij = distances) %>%
            filter(.data[['dij']] >= 0 & .data[['dij']] <= similarity.radii)
          if (nrow(info.neigh) == 0) {
            df <- df %>% mutate(flag_similarity = 1)
          } else {
            df.neigh <- df.neigh %>%
              filter(.data[[rid.column]] %in% info.neigh[[rid.column]])
            #
            mat.neigh <- tidyr::pivot_wider(df.neigh %>% select(all_of(c(date.column, rid.column, rdata.column))),
                                            names_from = all_of(rid.column),
                                            values_from = all_of(rdata.column))
            mat.target.neigh <- left_join(df.sub %>% select(all_of(c(date.column, data.column))), mat.neigh, by = date.column)
            mat.neigh <- mat.target.neigh %>% select(-all_of(c(date.column, data.column)))
            mat.target <- replicate(ncol(mat.neigh), mat.target.neigh[[data.column]])
            mat.similarity <- (mat.target - mat.neigh)^2
            similarity.stat <- apply(mat.similarity, 2, function(x) {
              zoo::rollapply(x, similarity.duration, FUN = function(y) sqrt(sum(y, na.rm = TRUE)), fill = NA_real_, align = 'right')
            })
            completeness.stat.neigh <- apply(mat.neigh, 2, function(x) {
              zoo::rollapply(x, similarity.duration, FUN = function(y) sum(!is.na(y)), fill = NA_real_, align = 'right')
            })
            completeness.stat.target <- zoo::rollapply(mat.target.neigh[[data.column]], similarity.duration,
                                                       FUN = function(y) sum(!is.na(y)), fill = NA_real_, align = 'right')
            similarity.stat <- ifelse(completeness.stat.neigh < compl.similarity.threshold, NA_real_, similarity.stat)
            completeness.stat.neigh <- ifelse(completeness.stat.neigh < compl.similarity.threshold, 0, 1)
            completeness.flag.target <- ifelse(completeness.stat.target < compl.similarity.threshold, 0, 1)
            completeness.flag.neigh <- apply(completeness.stat.neigh, 1, function(x) {
              ifelse(sum(is.na(x)) == length(x), NA_real_, ifelse(sum(!is.na(x)) >= 1, 1, 0))
            })
            #
            df.expected.similarity <- mat.target.neigh %>%
              mutate(Season = lubridate::quarter(.data[[date.column]], fiscal_start = 12)) %>%
              select(-all_of(data.column)) %>%
              tidyr::pivot_longer(cols = -all_of(c(date.column, "Season")),
                                  names_to = rid.column,
                                  values_to = 'value') %>%
              select(-.data[['value']])
            df.expected.similarity <- left_join(df.expected.similarity, info.neigh %>% select(all_of(c(rid.column, 'dij'))), by = rid.column) %>%
              mutate(expected_similarity = expected_similarity(season = .data[['Season']], distance = .data[['dij']] / 1000, relax.factor = 0.25))
            mat.expected.similarity <- pivot_wider(df.expected.similarity %>% select(all_of(c(date.column, rid.column, 'expected_similarity'))),
                                                   names_from = all_of(rid.column), values_from = expected_similarity) %>%
              select(-all_of(date.column))
            expected.similarity.stat <- apply(mat.expected.similarity, 2, function(x) {
              zoo::rollapply(x, similarity.duration, FUN = function(y) min(Mode_fun(y)), fill = NA_real_, align = 'right')
            })
            #
            mat.comparison <- similarity.stat < expected.similarity.stat
            similarity.flag <- apply(mat.comparison, 1, function(x) {
              if (all(is.na(x))) NA_real_ else if (sum(x, na.rm = TRUE) >= 1) 1 else 0
            })
            mat.comparison[is.na(mat.comparison) | mat.comparison == FALSE] <- NA_real_
            flag.similarity <- apply(mat.comparison, 1, function(x) {
              if (sum(x, na.rm = TRUE) >= 1) 3 else 2
            })
            flag.similarity <- ifelse(completeness.flag.target == 0, 0,
                                      ifelse(completeness.flag.neigh == 0, 1, flag.similarity))
            df.sub <- df.sub %>% mutate(flag_similarity = flag.similarity)
            df.out <- df %>% filter(.data[[date.column]] < start.nona | .data[[date.column]] > end.nona)
            df <- bind_rows(df.out, df.sub) %>%
              arrange(.data[[date.column]]) %>%
              mutate(flag_similarity = ifelse(is.na(.data[[data.column]]), NA_real_, .data[['flag_similarity']]))
          }
        }
      }
    }
    #
    if (new) {
      df <- df %>%
        mutate(!!paste0(param, "_similarity") := ifelse(.data[['flag_similarity']] != 3 | is.na(.data[['flag_similarity']]), NA_real_, .data[[data.column]]))
    }
    return(as_tibble(df))
  }
  #
  if (!is.null(id.column)) {
    data <- data %>% group_by(across(all_of(id.column))) %>%
      group_map(~ similarity_fun(.x), .keep = TRUE) %>%
      bind_rows()
  } else {
    data <- similarity_fun(data)
  }
  return(as_tibble(data))
}
