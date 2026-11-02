#' Correction Process
#'
#' Corrects the measurements at a location using data from nearby official monitoring stations and, optionally, auxiliary information (e.g., meteorological parameters), and assigns a quality flag indicating the success of the correction. Returns the corrected values together with their prediction intervals.
#'
#'
#' @details
#' Let \eqn{X_{it}} denote the measurements at a location \eqn{i} and \eqn{Y_{it}} the aggregated reference data from the nearby official monitoring stations.
#'
#' The correction model can be trained on a different version of the measurements than the one it is applied to (`data.flag.column`):
#' \itemize{
#'   \item if `data.flag.column = NULL` (default), the model is trained on and applied to the measurements in `data.column`, so that \eqn{X_{it}^{\ast} = X_{it}};
#'   \item otherwise, the model is trained on the flagged measurements \eqn{X_{it}^{\ast}} in `data.flag.column` (e.g., the measurements that passed the quality control, with `NA` elsewhere) and applied to the raw measurements \eqn{X_{it}} in `data.column`. In this way, the model is developed only with reliable measurements, while all the raw measurements are corrected.
#' }
#'
#' The correction can be applied to:
#' \itemize{
#'   \item the entire dataset (`historical.corr = FALSE`), where \eqn{\mathcal{T}} is the whole time range, or
#'   \item consecutive rolling windows (`historical.corr = TRUE`). One window is evaluated per day \eqn{n}, ending at 23:00 of that day and extending backward for `ndays.corr` days:
#'   \deqn{
#'   \mathcal{T}_n = [t_n - (T - 1)\Delta t, t_n],
#'   }{T_n = [t_n - (T - 1) dt, t_n]}
#'   where \eqn{t_n} is 23:00 of day \eqn{n}, \eqn{\Delta t} is the temporal resolution (1 hour) and \eqn{T = n_{\text{days}} \cdot 24} (`ndays.corr`) is the duration of the window in hours.
#' }
#'
#' The correction model is trained and evaluated either:
#' \itemize{
#'   \item with a split (`split = TRUE`): \eqn{\mathcal{T}_n} is divided into a training set \eqn{\mathcal{T}_n^{\ast}}, which contains all days except the last one and is used for the model fitting, and a test set \eqn{\mathcal{T}_n \setminus \mathcal{T}_n^{\ast}}, the last day, to which the model is applied, or
#'   \item without a split (`split = FALSE`): the whole window is used, so that \eqn{\mathcal{T}_n^{\ast} = \mathcal{T}_n}.
#' }
#'
#' The correction is performed in four steps:
#'
#' \strong{Step 1: Extraction of the nearby reference and meteorological stations}
#'
#' For each location \eqn{i}, the reference and meteorological stations within the maximum radius `m.radii` are identified with \code{\link{find_neighbors}}. The reference stations are then restricted to those within a season-dependent radius of influence \eqn{\rho_{\text{corr}}}{rho_corr}, obtained from a look-up table (see \code{\link{lut_distance_fun}}). The meteorological stations are used within `m.radii`.
#'
#' \strong{Step 2: Weighted reference and meteorological data}
#'
#' The reference data are spatially aggregated with an Inverse Distance Weighting scheme within \eqn{\rho_{\text{corr}}}{rho_corr}, and the auxiliary parameters within `m.radii`. See \code{\link{modified_idw_function}}.
#'
#' \strong{Step 3: Extraction of low-variability periods}
#'
#' Periods of low variability in the reference data are selected from the training set, for robust model development, with \code{\link{min_variability_function}} (`nmin`). Auxiliary parameters with insufficient data in the training set are not used.
#'
#' \strong{Step 4: Development of the correction model}
#'
#' A multiple linear regression model is trained on the selected data. The reference values are the dependent variable, while the measurements \eqn{X_{it}^{\ast}} at location \eqn{i} and the auxiliary parameters are the independent variables. The model is then applied to the measurements \eqn{X_{it}} of the test set and returns the corrected values \eqn{\hat{Y}_{it}} and their prediction intervals at the level \eqn{100(1-\alpha)\%} (`conf.level`). See \code{\link{correction_function}}.
#'
#' The performance of the model is assessed in the training period, by comparing the measurements used in the training and their corrected values with the reference values through the Root Mean Square Deviation (RMSD):
#' \deqn{
#' RMSD_{\text{raw}} = \sqrt{\frac{1}{N_{t^{\ast}}} \sum_{t^{\ast}} \left( Y_{it^{\ast}} - X_{it^{\ast}}^{\ast} \right)^2},
#' }{RMSD_raw = sqrt(mean((Y - X*)^2))}
#' \deqn{
#' RMSD_{\text{cor}} = \sqrt{\frac{1}{N_{t^{\ast}}} \sum_{t^{\ast}} \left( Y_{it^{\ast}} - \hat{Y}_{it^{\ast}} \right)^2},
#' }{RMSD_cor = sqrt(mean((Y - Yhat)^2))}
#' where \eqn{t^{\ast}} are the timestamps of the training set \eqn{\mathcal{T}_n^{\ast}} used in the model, and \eqn{N_{t^{\ast}}} their number.
#'
#' The quality flag \eqn{f_{i}^{\text{correction}}} is then defined as:
#' \deqn{
#' f_{i}^{\text{correction}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: } N_{\mathcal{T}_n^{\ast}} < n_{\text{min}}, \\
#' 1 & \text{Isolated site: } N_{\mathcal{J}}(\rho_{\text{corr}}) < 1, \\
#' 2 & \text{No improvement: } RMSD_{\text{cor}} > RMSD_{\text{raw}}, \\
#' 3 & \text{Improvement: } RMSD_{\text{cor}} \leq RMSD_{\text{raw}},
#' \end{cases}
#' }{f = 0 (insufficient data), 1 (isolated site), 2 (no improvement), 3 (improvement)}
#' \itemize{
#'   \item \eqn{\mathcal{J}}: set of the reference stations within \eqn{\rho_{\text{corr}}}{rho_corr}.
#'   \item \eqn{N_{\mathcal{J}}(\rho_{\text{corr}}) = \left| \{ j : d_{ij} \leq \rho_{\text{corr}} \text{ and } N_{j\mathcal{T}_n^{\ast}} \geq n_{\text{min}} \} \right|}{N_J = number of j with d_ij <= rho_corr and N_j >= n_min}: number of reference stations within \eqn{\rho_{\text{corr}}}{rho_corr} with at least \eqn{n_{\text{min}}} valid measurements in the training period.
#'   \item \eqn{N_{\mathcal{T}_n^{\ast}}}: number of valid measurements \eqn{X_{it}^{\ast}} of location \eqn{i} available for the model fitting (after the selection of the low-variability periods).
#'   \item \eqn{N_{j\mathcal{T}_n^{\ast}}}: number of valid measurements of the reference station \eqn{j} in the training period.
#'   \item \eqn{n_{\text{min}}}{n_min} (`compl.correction.threshold`): minimum number of non-missing measurements required to run the correction.
#' }
#'
#' When the correction is applied in consecutive windows (`historical.corr = TRUE`), the flag obtained for the window of day \eqn{n} is assigned to all timestamps within that day, \eqn{t \in n}:
#' \deqn{
#' f_{it}^{\text{correction}}(n) = f_{i}^{\text{correction}}.
#' }{f_it(n) = f_i}
#' Otherwise, the same flag is assigned to all timestamps:
#' \deqn{
#' f_{it}^{\text{correction}} = f_{i}^{\text{correction}}.
#' }{f_it = f_i}
#'
#' @param data A data frame containing the dataset to be corrected.
#' @param date.column character. Name of the date/time column in `data`, `rdata` and `mdata`.
#' @param data.column character. Name of the column in `data` to correct (e.g., the raw measurements).
#' @param data.flag.column character. Name of the column in `data` with the measurements used to train the correction model (e.g., the measurements that passed the quality control, with `NA` elsewhere). If `NULL`, the model is trained on `data.column`. Default is `NULL`. See Details.
#' @param id.column character. Name of the column identifying the measurement location ID in `data`.
#' @param rdata A data frame containing the dataset of the reference stations.
#' @param rdata.column character. Name of the column in `rdata` used in the correction.
#' @param rid.column character. Name of the column identifying the measurement location ID in `rdata`.
#' @param mdata A data frame containing the auxiliary information (e.g., meteorological parameters). Optional. Default is `NULL`.
#' @param mdata.column character. Vector of the names of the auxiliary parameter columns in `mdata`.
#' @param mid.column character. Name of the column identifying the measurement location ID in `mdata`.
#' @param x.column character. Name of the longitude column in `data`, `rdata` and `mdata`.
#' @param y.column character. Name of the latitude column in `data`, `rdata` and `mdata`.
#' @param m.radii numeric. Maximum radius of influence, in km, for extracting the reference and meteorological stations. The reference stations are then restricted to the season-dependent radius \eqn{\rho_{\text{corr}}}{rho_corr}, while the meteorological stations are used within `m.radii`. Default is `30`.
#' @param nmin numeric. Number of selected low-variability periods. See \code{\link{min_variability_function}}. Default is `5`.
#' @param compl.correction.threshold numeric. Minimum number of non-missing measurements required to run the correction. Default is `90`.
#' @param conf.level numeric. Confidence level of the prediction intervals of the corrected values. Default is `0.95`. See \code{\link{correction_function}}
#' @param robust logical. If `TRUE`, robust regression (`MASS::rlm`) is used; otherwise, ordinary least squares (`lm`). Default is `TRUE`. See \code{\link{correction_function}}
#' @param historical.corr logical. If `TRUE`, the correction is applied in consecutive daily rolling windows; otherwise, it is applied to the entire `data`. Default is `TRUE`.
#' @param ndays.corr numeric. Length of the rolling window, in days, if `historical.corr = TRUE`. Default is `30`.
#' @param param character. Base name of the new columns (e.g. `PM2.5` gives `PM2.5_cor`, `PM2.5_cor_lower` and `PM2.5_cor_upper`). Default is `"PM2.5"`.
#' @param split logical. If `TRUE`, the data are divided into training and test sets for the development and the application of the model. Default is `TRUE`. See Details.
#'
#' @return The original dataset with the following columns added:
#' \itemize{
#'   \item \code{<param>_cor}: numeric. Corrected values of `data.column`.
#'   \item \code{<param>_cor_lower}: numeric. Lower bound of the prediction interval of the corrected values (e.g., at the 95\% level, the 2.5\% percentile).
#'   \item \code{<param>_cor_upper}: numeric. Upper bound of the prediction interval of the corrected values (e.g., at the 95\% level, the 97.5\% percentile).
#'   \item \code{rdata.column}: numeric. Reference observations at the sensor location, interpolated from the nearby reference stations using IDW, used to derive and evaluate the correction.
#'   \item \code{flag_correction}: integer. Quality flags for the correction process, \code{0} = insufficient data, \code{1} = isolated site, \code{2} = no improvement of the initial data, \code{3} = improvement of the initial data. It is `NA` where the value to correct is missing.
#' }
#'
#' @seealso \code{\link{find_neighbors}}, \code{\link{modified_idw_function}}, \code{\link{min_variability_function}}, and \code{\link{correction_function}}
#'
#' @export
correction_process <- function(data,
                               date.column,
                               data.column,
                               data.flag.column = NULL,
                               id.column,
                               rdata,
                               rdata.column,
                               rid.column,
                               mdata,
                               mdata.column,
                               mid.column,
                               x.column,
                               y.column,
                               m.radii = 30,
                               nmin = 5,
                               compl.correction.threshold = 90,
                               conf.level = 0.95,
                               robust = TRUE,
                               historical.corr = TRUE,
                               ndays.corr = 30,
                               param = 'PM2.5',
                               split = TRUE) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(data), date.column %in% colnames(rdata))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  if (!is.null(data.flag.column)) {
    stopifnot(is.character(data.flag.column), data.flag.column %in% colnames(data))
  }
  # column used to train the correction: the flagged data (if given), otherwise the raw data
  train.column <- if (is.null(data.flag.column)) data.column else data.flag.column
  stopifnot(is.character(param), length(param) == 1)
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(id.column), id.column %in% colnames(data))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(x.column),x.column %in% colnames(data), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column),y.column %in% colnames(data), y.column %in% colnames(rdata))
  stopifnot(is.numeric(compl.correction.threshold), !is.null(compl.correction.threshold))
  stopifnot(is.numeric(ndays.corr),!is.null(ndays.corr))
  stopifnot(is.logical(robust))
  stopifnot(is.logical(historical.corr), is.logical(split))
  conf.level <- ifelse(is.null(conf.level) | is_missing(conf.level), 0.95, conf.level)
  nmin <- ifelse(is.null(nmin), 24, nmin)
  #
  if (id.column == rid.column) {
    rdata <- rdata %>%
      rename(!!paste0('r_', rid.column) := rid.column)
    rid.column <- paste0('r_',rid.column)
  }
  #
  data %>%
    group_by(across(all_of(id.column))) %>%
    group_map(~ {
      id.target <- .y[[1]]
      df <- .x
      # all proximal reference stations (up to m.radii)
      df.target.all.ref <- find_neighbors(data=df,
                                          rdata=rdata,
                                          date.column = date.column,
                                          data.column = data.column,
                                          rdata.column = rdata.column,
                                          id.column = id.column,
                                          rid.column = rid.column,
                                          x.column = x.column,
                                          y.column = y.column,
                                          maxradii = m.radii,
                                          maxneigh = NULL,
                                          compl.threshold = compl.correction.threshold)
      df.nearest.all.ref <- df.target.all.ref$neighbors
      df.target.all.ref <- df.target.all.ref$data
      # all proximal meteorological stations (up to m.radii)
      df.target.all.met <- find_neighbors(data=df,
                                          rdata=mdata,
                                          date.column = date.column,
                                          data.column = data.column,
                                          rdata.column = mdata.column,
                                          id.column = id.column,
                                          rid.column = mid.column,
                                          x.column = x.column,
                                          y.column = y.column,
                                          maxradii = m.radii,
                                          maxneigh = NULL,
                                          compl.threshold = compl.correction.threshold)
      df.nearest.all.met <- df.target.all.met$neighbors
      df.target.all.met <- df.target.all.met$data
      #
      correction_process_fun <- function(df.target, split = split) {
        df.target.compl <- df.target %>%
          tidyr::drop_na(all_of(train.column)) %>%
          dplyr::summarise(n_hours = sum(!is.na(.data[[train.column]])))
        correction.radii <- lut_distance_fun(season=min(Mode_fun(quarter(df.target[[date.column]], fiscal_start = 12))))
        #
        if (split) {
          df.target <- df.target %>%
            mutate(day = as_date(.data[[date.column]])) %>%
            arrange(.data[[date.column]])
          df.target.train <- df.target %>%
            filter(day < max(day)) %>%
            select(-day)
          df.target.test <- df.target %>%
            filter(day == max(day)) %>%
            select(-day)
          df.target <- df.target %>%
            select(-day)
        } else {
          df.target.train <- df.target.test <- df.target
        }
        # the correction is trained on the flagged data and applied to the raw data:
        # in the training data set, the raw values are replaced by the flagged values
        if (!is.null(data.flag.column)) {
          df.target.train[[data.column]] <- df.target.train[[data.flag.column]]
        }
        #
        if (df.target.compl[['n_hours']] < compl.correction.threshold) {
          df.target.test <- df.target.test %>%
            mutate(flag_correction = 0,
                   !!(paste0(param,'_cor')) := NA_real_,
                   !!(paste0(param,'_cor_lower')) := NA_real_,
                   !!(paste0(param,'_cor_upper')) := NA_real_,
                   !!rdata.column := NA_real_)
        } else {
          # proximal reference station up to correction.radii
          df.nearest.ref <- df.nearest.all.ref %>%
            filter(dij <= correction.radii)
          df.target.ref <- df.target.all.ref %>%
            filter(.data[[rid.column]] %in% df.nearest.ref[[rid.column]])
          #
          if (is.null(df.nearest.ref)) {
            df.target.test <- df.target.test %>%
              mutate(flag_correction = 1,
                     !!(paste0(param,'_cor')) := NA_real_,
                     !!(paste0(param,'_cor_lower')) := NA_real_,
                     !!(paste0(param,'_cor_upper')) := NA_real_,
                     !!rdata.column := NA_real_)
          }
          df.weighted.ref <- modified_idw_function(info = df.nearest.ref,
                                                   data = df.target.ref  %>% dplyr::select(-any_of(c(id.column, data.column, data.flag.column))),
                                                   radii = correction.radii,
                                                   date.column = date.column,
                                                   id.column = rid.column,
                                                   data.column = rdata.column)
          df.target.ref <- left_join(df.target, df.weighted.ref, by=date.column)
          #
          if (!is.null(mdata)) {
            stopifnot(is.character(date.column), date.column %in% colnames(mdata))
            stopifnot(is.character(mdata.column), mdata.column %in% colnames(mdata))
            stopifnot(is.character(x.column), x.column %in% colnames(mdata))
            stopifnot(is.character(y.column), y.column %in% colnames(mdata))
            stopifnot(is.character(mid.column), mid.column %in% colnames(mdata))
            m.radii <- ifelse(is.null(m.radii), 30, m.radii)
            #
            if (id.column == mid.column | rid.column == mid.column) {
              mdata <- mdata %>%
                rename(!!paste0('m_', mid.column) := mid.column)
              mid.column <- paste0('m_',mid.column)
            }
            # proximal meteorological stations up to m.radii
            df.nearest.met <- df.nearest.all.met %>%
              filter(dij <= m.radii)
            df.target.met <- df.target.all.met %>%
              filter(.data[[mid.column]] %in% df.nearest.met[[mid.column]])
            #
            if (!is.null(df.nearest.met)) {
              df.target.met <- df.target.met %>%
                group_by(.data[[mid.column]]) %>%
                mutate(across(all_of(mdata.column), ~ na.approx(., na.rm=FALSE))) %>%
                ungroup()
              #
              df.weighted.met <- modified_idw_function(info = df.nearest.met,
                                                       data = df.target.met  %>% dplyr::select(-any_of(c(id.column, data.column, data.flag.column))),
                                                       radii = m.radii,
                                                       date.column = date.column,
                                                       id.column = mid.column,
                                                       data.column = mdata.column)
              df.target.ref <- left_join(df.target.ref, df.weighted.met, by=date.column) %>% as_tibble()
            }
          }
          #
          # train-test datasets
          df.target.train <- left_join(df.target.train, df.target.ref %>% select(any_of(c(date.column, rdata.column, mdata.column))), by=date.column)
          df.target.test <- left_join(df.target.test, df.target.ref %>% select(any_of(c(date.column, rdata.column, mdata.column))), by=date.column)
          #
          # check the availability of meteorological data and remove meteorological variables with non-adequate data (in the training data set)
          meteo.compl <- df.target.train %>%
            drop_na(any_of(c(data.column, rdata.column))) %>%
            summarise(across(any_of(mdata.column), ~ sum(!is.na(.)))) %>%
            select(where(~ . > compl.correction.threshold))
          mdata.keep <- names(meteo.compl)
          #
          #  training - calibration data
          # calculation of the minimum reference values
          df.target.train <- min_variability_function(data=df.target.train, date.column = date.column, data.column = rdata.column, nmin = nmin, crop = TRUE)
          df.target.train <- df.target.train %>%
            tidyr::drop_na(all_of(c(data.column, rdata.column)))
          #
          if (nrow(df.target.train) < compl.correction.threshold) {
            df.target.test <- df.target.test %>%
              dplyr::mutate(flag_correction = 0,
                            !!(paste0(param,'_cor')) := NA_real_,
                            !!(paste0(param,'_cor_lower')) := NA_real_,
                            !!(paste0(param,'_cor_upper')) := NA_real_,
                            !!rdata.column := NA_real_)
          } else {
            # correction process
            df.target.test <- df.target.test %>%
              mutate(across(any_of(mdata.keep), ~ if (sum(!is.na(.)) > 5) {na.approx(., na.rm=FALSE)} else {rep(NA_real_, length(.))}))
            mdata.keep <- df.target.test %>%
              drop_na(any_of(mdata.keep)) %>%
              summarise(across(any_of(mdata.keep), ~ sum(!is.na(.)))) %>%
              select(where(~ . > 5)) %>%
              names()
            #
            calibration.model <- correction_function(df = df.target.train,
                                                     y = rdata.column,
                                                     x = data.column,
                                                     covariates = if (length(mdata.keep) == 0) NULL else mdata.keep,
                                                     conf.level = conf.level,
                                                     newdata = df.target.test,
                                                     param = param,
                                                     robust = robust)
            
            df.target.train <- calibration.model$target
            df.target.test <- calibration.model$new
            df.calib.stats <- tibble(RMSDraw = rmsd(df.target.train[[rdata.column]], df.target.train[[data.column]]),
                                     RMSDcor = rmsd(df.target.train[[rdata.column]], df.target.train[[paste0(param, '_cor')]]))
            # calibration  flag
            flag_correction <- ifelse(sum(is.na(df.target.test[[data.column]])) == nrow(df.target.test), NA_real_,
                                      ifelse(is.na(df.calib.stats[['RMSDcor']]) & !is.na(df.calib.stats[['RMSDraw']]), 2,
                                             ifelse(df.calib.stats[['RMSDcor']] > df.calib.stats[['RMSDraw']], 2, 3)))
            df.target.test <- df.target.test %>%
              dplyr::mutate(flag_correction = flag_correction)
          }
        }
        #
        df.target.test <- df.target.test %>%
          dplyr::select(all_of(c(date.column, paste0(param,'_cor'), paste0(param,'_cor_lower'),paste0(param,'_cor_upper'), rdata.column, 'flag_correction')))
        #
        return(as_tibble(df.target.test))
      }
      #
      # check for (not) historical data
      if (historical.corr) {
        # run in temporal windows
        days.intervals <- ymd_hms(paste(unique(as_date(df[[date.column]])), "23:00:00"), tz= unique(tz(df[[date.column]])))
        days <- as_date(floor_date(days.intervals, unit = 'day'))
        df.correction.list <- runner::runner(
          x = df,
          f = function(x) correction_process_fun(x, split=split),
          k = paste(ndays.corr, "days"),
          idx = df[[date.column]],
          at = days.intervals,
          simplify = FALSE
        )
        # update for df
        df.correction <- bind_rows(df.correction.list)
        df <- left_join(df, df.correction, by = date.column)
        df <- df %>%
          mutate(flag_correction := ifelse(is.na(.data[[data.column]]), NA_real_, flag_correction))
      } else {
        df.correction <- correction_process_fun(df, split = split)
        df <- left_join(df, df.correction, by = date.column)
        df <- df %>%
          mutate(flag_correction := ifelse(is.na(.data[[data.column]]), NA_real_, flag_correction))
      }
      return(as_tibble(df))
    }, .keep = TRUE) %>%
    bind_rows()
}