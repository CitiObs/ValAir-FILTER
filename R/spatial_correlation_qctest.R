#' Spatial Correlation Test
#'
#' Checks whether the measurements at a location are correlated with those of neighboring locations, using Pearson's R correlation coefficient, and assigns a quality flag. Optionally creates a column of valid values.
#'
#'
#' @details
#' Let \eqn{X_{it}} denote the measurements at a location \eqn{i}. The test can be applied to:
#' \itemize{
#'   \item the entire dataset (`historical = FALSE`), where \eqn{\mathcal{T}} is the whole time range, or
#'   \item consecutive rolling windows (`historical = TRUE`). One window is evaluated per day \eqn{n}, ending at 23:00 of that day and extending backward for `ndays` days:
#'   \deqn{
#'   \mathcal{T} = [t_n - (T - 1)\Delta t, t_n],
#'   }
#'   where \eqn{t_n} is 23:00 of day \eqn{n}, \eqn{\Delta t} is the temporal resolution (1 hour) and \eqn{T = n_{\text{days}} \cdot 24} (`ndays`) is the duration of the window in hours.
#' }
#'
#' Missing values are excluded from the calculation. Let \eqn{\mathcal{T}^{\ast} \subseteq \mathcal{T}} be the set of timestamps in the window with available (non-missing) values for both \eqn{X_{it}} and \eqn{X_{jt}}, and \eqn{N_{\mathcal{T}^{\ast}}} their number. For each location \eqn{i}, the pairwise Pearson's R correlation coefficient \eqn{R_{ij}} is calculated between \eqn{i} and each neighboring location \eqn{j}, together with their separation distance \eqn{d_{ij}}:
#' \deqn{
#' R_{ij}(\mathcal{T}) = \frac{\sum\limits_{t \in \mathcal{T}^{\ast}} (X_{it} - \bar{X}_i^{\ast})(X_{jt} - \bar{X}_j^{\ast})}{\sqrt{\sum\limits_{t \in \mathcal{T}^{\ast}} (X_{it} - \bar{X}_i^{\ast})^2} \sqrt{\sum\limits_{t \in \mathcal{T}^{\ast}} (X_{jt} - \bar{X}_j^{\ast})^2}},
#' }{R_ij = Pearson correlation between X_it and X_jt over T*}
#' \itemize{
#'   \item \eqn{\bar{X}_i^{\ast}} and \eqn{\bar{X}_j^{\ast}}: temporal averages of \eqn{X_{it}} and \eqn{X_{jt}} within \eqn{\mathcal{T}^{\ast}}:
#'   \deqn{
#'   \bar{X}_i^{\ast} = \frac{1}{N_{\mathcal{T}^{\ast}}} \sum_{t \in \mathcal{T}^{\ast}} X_{it}, \quad \bar{X}_j^{\ast} = \frac{1}{N_{\mathcal{T}^{\ast}}} \sum_{t \in \mathcal{T}^{\ast}} X_{jt}
#'   }{mean of X_it and X_jt over T*}
#' }
#'
#' The calculated \eqn{R_{ij}(\mathcal{T})} is compared with an expected correlation \eqn{R_{\text{exp}}(d_{ij}, s)}, derived from a look-up table based on the separation distance \eqn{d_{ij}} and the season \eqn{s} of \eqn{t \in \mathcal{T}}. See \code{\link{expected_correlation}} for more information.
#'
#' For each window \eqn{\mathcal{T}}, the number of neighbors \eqn{j} within a radius of influence \eqn{\rho} that satisfy the following conditions is derived:
#' \itemize{
#'   \item Adequate pairwise measurements:
#'   \deqn{
#'   N_{j}(\rho) = \left| \{ j : N_{\mathcal{T}^{\ast}} \geq n_{\text{min}} \text{ and } d_{ij} \leq \rho \} \right|
#'   }{N_j = number of j with N_T* >= n_min and d_ij <= rho}
#'   \item Correlation above the expected threshold, \eqn{R_{ij}(\mathcal{T}) > R_{\text{exp}}(d_{ij}, s)}:
#'   \deqn{
#'   N_{R_{ij} > R_{\text{exp}}}(\rho) = \left| \{ j : R_{ij} > R_{\text{exp}} \text{ and } d_{ij} \leq \rho \} \right|
#'   }{N_corr = number of j with R_ij > R_exp and d_ij <= rho}
#' }
#'
#' The quality flag \eqn{f_{i}^{\text{correlation}}} for the location \eqn{i} is then defined as:
#' \deqn{
#' f_{i}^{\text{correlation}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: } N_{i} < n_{\text{min}}, \\
#' 1 & \text{Isolated site: } N_{j}(\rho_{\text{local}}) < 1 \text{ and } N_{j}(\rho_{\text{regional}}) < n_{\text{min,neigh}}, \\
#' 2 & \text{No correlation: } N_{R_{ij} > R_{\text{exp}}}(\rho_{\text{local}}) < 1 \text{ and } N_{R_{ij} > R_{\text{exp}}}(\rho_{\text{regional}}) < n_{\text{min,neigh}}, \\
#' 3 & \text{Correlation: } N_{R_{ij} > R_{\text{exp}}}(\rho_{\text{local}}) \geq 1 \text{ or } N_{R_{ij} > R_{\text{exp}}}(\rho_{\text{regional}}) \geq n_{\text{min,neigh}},
#' \end{cases}
#' }{f = 0 (insufficient data), 1 (isolated site), 2 (no correlation), 3 (correlation)}
#' \itemize{
#'   \item \eqn{N_{i}}: number of available (non-missing) values of location \eqn{i} within the window.
#'   \item \eqn{n_{\text{min}}}{n_min} (`compl.corr.threshold`): minimum number of non-missing values required to evaluate the test, given either as a number of measurements or as a fraction \eqn{p} of the window duration \eqn{T}, in which case \eqn{n_{\text{min}} = p T}{n_min = p * T}.
#'   \item \eqn{\rho_{\text{local}}}{rho_local}, \eqn{\rho_{\text{regional}}}{rho_regional} (`corr.radii`): inner and outer radii of influence, in km.
#'   \item \eqn{n_{\text{min,neigh}}}{n_min,neigh} (`minneigh`): minimum number of neighbors required within \eqn{\rho_{\text{regional}}}{rho_regional}.
#' }
#'
#' When the test is applied in consecutive windows (`historical = TRUE`), the flag obtained for the window of day \eqn{n} is assigned to all timestamps within that day, \eqn{t \in n}:
#' \deqn{
#' f_{it}^{\text{correlation}}(n) = f_{i}^{\text{correlation}}.
#' }{f_it(n) = f_i}
#' Otherwise, the same flag is assigned to all timestamps:
#' \deqn{
#' f_{it}^{\text{correlation}} = f_{i}^{\text{correlation}}.
#' }{f_it = f_i}
#'
#' @param data A data frame containing the dataset to be tested.
#' @param date.column character. Name of the date/time column in `data` and `rdata`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location ID in `data`.
#' @param id.target character. ID of a measurement location. If `NULL`, the correlation is calculated for all IDs in `id.column`. Default is `NULL`.
#' @param rdata A data frame containing the dataset of additional locations (e.g., reference stations) used as neighbors.
#' @param rdata.column character. Name of the column in `rdata` used in the test.
#' @param rid.column character. Name of the column identifying the measurement location ID in `rdata`.
#' @param x.column character. Name of the longitude column in `data` and `rdata`.
#' @param y.column character. Name of the latitude column in `data` and `rdata`.
#' @param compl.corr.threshold numeric. Minimum data completeness required within the window to run the test, expressed either as a number of measurements or as a fraction of the window duration. Default is `90`.
#' @param corr.radii numeric. Vector of two values with the inner (local) and the outer (regional) radius of influence, in km. Default is `c(3, 30)`.
#' @param minneigh numeric. Minimum number of neighboring locations to test against. Default is `2`.
#' @param historical logical. If `TRUE`, the test is applied in consecutive daily rolling windows; otherwise, it is applied to the entire `data`. Default is `TRUE`.
#' @param ndays numeric. Length of the rolling window, in days, if `historical = TRUE`. Default is `30`.
#' @param new logical. If TRUE, adds a column `<param>_correlation` containing the valid values (flag `3`); other values are replaced with `NA`. Default is `FALSE`.
#' @param param character. Base name of the new column (e.g. `PM2.5` gives `PM2.5_correlation`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_correlation}: integer. Quality flags based on the statistical test, \code{0} = insufficient data, \code{1} = isolated site, \code{2} = no correlation with neighbors, \code{3} = correlation with neighbors. It is `NA` where the tested value is missing.
#'   \item \code{<param>_correlation}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
spatial_correlation_qctest <- function(data,
                                       date.column,
                                       data.column,
                                       id.column,
                                       id.target = NULL,
                                       rdata = NULL,
                                       rdata.column = NULL,
                                       rid.column = NULL,
                                       x.column,
                                       y.column,
                                       compl.corr.threshold = 90,
                                       corr.radii = c(3, 30),
                                       minneigh = 2,
                                       historical = TRUE,
                                       ndays = 30,
                                       new = FALSE,
                                       param = "PM2.5") {
  # Suppress warnings inside this function, restore the user's setting on exit
  old.options <- options(warn = -1)
  on.exit(options(old.options), add = TRUE)
  #
  # Argument checks
  stopifnot(all(sapply(c(date.column, data.column, id.column, x.column, y.column), is.character)))
  stopifnot(all(c(date.column, data.column, id.column, x.column, y.column) %in% colnames(data)))
  stopifnot(is.numeric(compl.corr.threshold) && compl.corr.threshold >= 0)
  stopifnot(is.numeric(corr.radii), !is.null(corr.radii))
  stopifnot(is.numeric(minneigh), !is.null(minneigh))
  stopifnot(is.logical(historical))
  #
  data <- data %>%
    mutate(!!id.column := as.character(.data[[id.column]]))
  # correlation radii to m
  corr.radii <- corr.radii * 1000
  #
  single_location_correlation_fun <- function(df, id.target) {
    correlation_process_fun <- function(df.target) {
      if (sum(!is.na(df.target[[data.column]])) < compl.corr.threshold) return(0)
      #
      # target location
      timestamps.range <- range(df.target[[date.column]])
      info.target <- df.target %>%
        filter(.data[[id.column]] == id.target) %>%
        distinct(.data[[id.column]], .data[[x.column]], .data[[y.column]])
      sp.target <- st_as_sf(info.target, coords = c(x.column, y.column), crs = 4326L)
      #
      # neighbor data and incorporate rdata if exists
      if (!is.null(rdata)) {
        # arguments for rdata
        stopifnot(all(sapply(list(rdata.column, date.column, rid.column), function(x) !is.null(x) && is.character(x))))
        stopifnot(all(c(rdata.column, date.column, rid.column) %in% colnames(rdata)))
        #
        rdata <- rdata %>%
          select(all_of(c(date.column, rid.column, rdata.column, x.column, y.column))) %>%
          rename(!!id.column := !!rid.column,
                 !!data.column := !!rdata.column) %>%
          mutate(!!id.column := paste0("ref_", .data[[id.column]])) %>%
          select(all_of(c(date.column, id.column, x.column, y.column, data.column)))
        df.neigh <- bind_rows(data %>% filter(.data[[id.column]] != id.target) , rdata) %>% as_tibble()
      } else {
        df.neigh <- data %>%
          filter(.data[[id.column]] != id.target)
      }
      df.neigh <- df.neigh %>%
        filter(.data[[date.column]] >= timestamps.range[1],
               .data[[date.column]] <= timestamps.range[2])
      #
      # omit neighbor data outside the temporal window
      df.neigh <- df.neigh %>%
        group_by(across(all_of(id.column))) %>%
        filter(sum(!is.na(.data[[data.column]])) >= compl.corr.threshold) %>%
        ungroup()
      #
      if (nrow(df.neigh) == 0) return(0)
      info.neigh <- df.neigh %>%
        distinct(.data[[id.column]], .data[[x.column]], .data[[y.column]])
      sp.neigh <- sf::st_as_sf(info.neigh, coords = c(x.column, y.column), crs = 4326L)
      # calculation of separation distance and selection of neighbors
      distances <- as.vector(round(st_distance(sp.target, sp.neigh), 1))
      info.neigh <- info.neigh %>%
        mutate(dij = distances) %>%
        filter(.data[["dij"]] < corr.radii[2])
      df.neigh <- df.neigh %>%
        filter(.data[[id.column]] %in% info.neigh[[id.column]])
      if (nrow(df.neigh) == 0) return(0)
      #
      # neighboring and target data in matrix format
      mat.neigh <- pivot_wider(df.neigh %>% select(all_of(c(date.column, id.column, data.column))),
                               names_from = all_of(id.column),
                               values_from = all_of(data.column))
      mat.neigh <- left_join(df.target %>% select(all_of(c(date.column, data.column))), mat.neigh, by = date.column) %>%
        select(-all_of(c(date.column, data.column)))
      #
      # Compute pairwise correlation between target and neighbors
      mat.target <- replicate(ncol(mat.neigh), df.target[[data.column]])
      mat.completeness <- apply(mat.target + mat.neigh, 2, function(x) sum(!is.na(x)))
      completeness.stat <- tibble(!!id.column := names(mat.completeness), completeness = as.vector(mat.completeness))
      pairwise.correlation <- mapply(cor, as.data.frame(mat.neigh), as.data.frame(mat.target), method='pearson', use='pairwise.complete.obs')
      correlation.stat <- tibble(!!id.column := names(pairwise.correlation), correlation = as.vector(pairwise.correlation))
      df.correlation <- left_join(correlation.stat, completeness.stat, by=id.column)
      #
      # Calculate expected correlation by season and distance
      df.expected.correlation <- df.neigh %>%
        mutate(Season = quarter(.data[[date.column]], fiscal_start = 12)) %>%
        select(all_of(c("Season", id.column))) %>%
        group_by(.data[[id.column]]) %>%
        filter(.data[["Season"]] == min(Mode_fun(.data[["Season"]]))) %>%
        dplyr::distinct() %>% ungroup()
      df.expected.correlation <- left_join(df.expected.correlation, info.neigh %>% select(all_of(c(id.column, "dij"))), by = id.column)
      df.expected.correlation <- df.expected.correlation %>%
        mutate(expected_correlation = expected_correlation(season = .data[["Season"]], distance = .data[["dij"]] / 1000, relax.factor = 0.25))
      #
      # Comparisons
      df.comparison <- left_join(df.correlation, df.expected.correlation %>% select(-all_of("Season")), by = id.column) %>%
        filter(.data[["completeness"]] > compl.corr.threshold)
      if (nrow(df.comparison) == 0) return(1)
      df.comparison <- df.comparison %>%
        filter(.data[["correlation"]] > .data[["expected_correlation"]])
      # Check neighborhood criteria for return code 3
      if (sum(df.comparison[["dij"]] <= corr.radii[1], na.rm = TRUE) >= 1 |
          sum(df.comparison[["dij"]] > corr.radii[1], na.rm = TRUE) >= minneigh) {
        return(3)
      }
      return(2)
    }
    #
    # check for (not) historical data
    if (historical) {
      stopifnot(is.numeric(ndays), ndays > 0, !is.null(ndays))
      # run in temporal windows
      days.intervals <- ymd_hms(paste(unique(as_date(df[[date.column]])), "23:00:00"))
      days <- as_date(floor_date(days.intervals, unit = 'day'))
      correlation.flags <- runner(x = df,
                                  f = function(x) correlation_process_fun(x),
                                  k = paste(ndays, "days"),
                                  idx = df[[date.column]],
                                  at = days.intervals,
                                  simplify = TRUE)
      # update for df
      daily.correlation.flag <- tibble(Date = days, flag_correlation = correlation.flags)
      df <- df %>%
        mutate(Date = as_date(.data[[date.column]])) %>%
        left_join(daily.correlation.flag, by = "Date") %>%
        select(-Date)
      df <- df %>%
        mutate(flag_correlation = ifelse(is.na(.data[[data.column]]), NA_real_, .data[["flag_correlation"]]))
    } else {
      correlation.flags <- correlation_process_fun(df)
      df <- df %>%
        mutate(flag_correlation = ifelse(is.na(.data[[data.column]]), NA_real_, correlation.flags))
    }
    if (new == TRUE) {
      df <- df %>%
        mutate(!!paste0(param, "_correlation") := ifelse(.data[["flag_correlation"]] != 3 | is.na(.data[["flag_correlation"]]), NA_real_, .data[[data.column]]))
    }
    return(as_tibble(df))
  }
  #
  # all vs. single location
  if (is.null(id.target)) {
    output <- data %>%
      group_by(across(all_of(id.column))) %>%
      group_map(
        ~ single_location_correlation_fun(df = .x, id.target = .y[[1]]), .keep = TRUE) %>%
      bind_rows()
  } else {
    ## Run only one target location
    id.target <- as.character(id.target)
    stopifnot(id.target %in% data[[id.column]])
    output <- single_location_correlation_fun(df = data %>% filter(.data[[id.column]] == id.target), id.target = id.target)
  }
  return(as_tibble(output))
}
