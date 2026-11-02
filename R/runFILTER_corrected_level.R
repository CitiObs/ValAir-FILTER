#' Quality Control based on FILTER (Corrected Processing Level)
#'
#' Applies the full quality control (QC) procedure of FILTER on the corrected data.
#'
#' @details
#' This function performs all quality control (QC) algorithms included in the FILTER approach on the corrected data. The correction model is trained using only the
#' measurements that passed the range validity, constant value and outlier detection (temporal and spatial) tests (see \code{\link{range_validity_qctest}},
#' \code{\link{constant_value_qctest}}, \code{\link{outlier_qctest}}, and \code{\link{spatial_outlier_qctest}}), i.e., with \code{flag_range = 1},
#' \code{flag_constant_value = 2} and \code{flag_outlier = 2}, and is then applied to all the raw measurements (see `data.flag.column` in \code{\link{correction_process}}).\cr\cr
#'
#' \strong{Raw Data} ---------  \strong{Range Validity} --> \strong{Constant Value} --> \strong{Outlier Detection} --> \strong{Correction process} --------- \strong{Corrected Data}\cr\cr
#'
#' Then, QC is applied on the corrected data using the five-step procedure of FILTER (see \code{\link{runFILTER_raw_level}}).\cr\cr
#'
#' \strong{Corrected Data} ---------  \strong{Range Validity} --> \strong{Constant Value} --> \strong{Outlier Detection} --> \strong{Spatial Correlation} --> \strong{Spatial Similarity} --------- \strong{QC Corrected Data}\cr\cr
#'
#' At the end of the process, a total quality flag is assigned to each corrected value.
#'
#' @param data A data frame containing the dataset to be tested.
#' @param date.column character. Name of the date/time column in `data`.
#' @param data.column character. Name of the column in `data` on which the test is performed.
#' @param id.column character. Name or ID of the measurement location in `data`.
#' @param upper.bound see in \code{\link{range_validity_qctest}}.
#' @param lower.bound see in \code{\link{range_validity_qctest}}.
#' @param persist.duration see in \code{\link{constant_value_qctest}}.
#' @param compl.persist.duration see in \code{\link{constant_value_qctest}}.
#' @param min.variation see in \code{\link{constant_value_qctest}}.
#' @param method see in \code{\link{constant_value_qctest}}.
#' @param outlier.duration see in \code{\link{outlier_qctest}}.
#' @param compl.outlier.duration see in \code{\link{outlier_qctest}}.
#' @param outlier.threshold see in \code{\link{outlier_qctest}}.
#' @param out.flag.column see in \code{\link{spatial_outlier_qctest}}.
#' @param outlier.radii see in \code{\link{spatial_outlier_qctest}}.
#' @param outlier.minneigh see in \code{\link{spatial_outlier_qctest}}.
#' @param rdata see in \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{correction_process}}.
#' @param rdata.column see in \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{correction_process}}.
#' @param rid.column see in \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{correction_process}}.
#' @param compl.corr.threshold see in \code{\link{spatial_correlation_qctest}}.
#' @param corr.radii see in \code{\link{spatial_correlation_qctest}}.
#' @param minneigh see in \code{\link{spatial_correlation_qctest}}.
#' @param historical see in \code{\link{spatial_correlation_qctest}}. Applies only to the Spatial Correlation test; for the correction, see `historical.corr`.
#' @param ndays see in \code{\link{spatial_correlation_qctest}}.
#' @param similarity.radii see in \code{\link{spatial_similarity_qctest}}.
#' @param similarity.duration see in \code{\link{spatial_similarity_qctest}}.
#' @param compl.similarity.threshold see in \code{\link{spatial_similarity_qctest}}.
#' @param mdata see in \code{\link{correction_process}}.
#' @param mdata.column see in \code{\link{correction_process}}.
#' @param mid.column see in \code{\link{correction_process}}.
#' @param m.radii see in \code{\link{correction_process}}.
#' @param historical.corr see in \code{\link{correction_process}}.
#' @param ndays.corr see in \code{\link{correction_process}}.
#' @param nmin see in \code{\link{correction_process}}.
#' @param compl.correction.threshold see in \code{\link{correction_process}}.
#' @param conf.level see in \code{\link{correction_function}} and \code{\link{correction_process}}.
#' @param robust see in \code{\link{correction_function}} and \code{\link{correction_process}}.
#' @param split see in \code{\link{correction_process}}.
#' @param qcflags see in \code{\link{quality_flag_function}}.
#' @param x.column see in \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, and \code{\link{spatial_similarity_qctest}}.
#' @param y.column see in \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, and \code{\link{spatial_similarity_qctest}}.
#' @param param see in \code{\link{range_validity_qctest}}, \code{\link{constant_value_qctest}}, \code{\link{outlier_qctest}}, \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{correction_process}}.
#' @param keep character. Which quality flags of the corrected data are kept in the output. If `"flags"` (default), the individual quality flags of each test are kept together with the total quality flag (`quality_flag`). If `"none"`, the individual quality flags are removed, so only the total quality flag is added. See \code{\link{runFILTER_raw_level}}.
#'
#' @return Returns the original data including the corrected values with their confidence intervals (see \code{\link{correction_process}}) and the total quality flag (`quality_flag`) of each corrected value. The flag has three levels: `2` (high quality), `1` (good quality) and `0` (other quality). If `keep = "flags"`, the individual quality flags of each test on the corrected data are also included.
#'
#' @seealso \code{\link{correction_process}} and \code{\link{runFILTER_raw_level}}.
#'
#' @export
runFILTER_corrected_level <- function(data,
                                      date.column,
                                      data.column,
                                      id.column,
                                      # range test arguments
                                      upper.bound = 1000,
                                      lower.bound = 0,
                                      # constant value test arguments
                                      persist.duration = 8,
                                      compl.persist.duration = 0.75,
                                      min.variation = 0.1,
                                      method = 'range',
                                      # outlier test arguments
                                      outlier.duration = 168,
                                      compl.outlier.duration = 0.25,
                                      outlier.threshold = 10,
                                      # spatial outlier arguments
                                      out.flag.column = 'flag_outlier',
                                      outlier.radii = c(3, 30),
                                      outlier.minneigh = 2,
                                      # spatial correlation arguments
                                      rdata,
                                      rdata.column,
                                      rid.column,
                                      compl.corr.threshold = 90,
                                      corr.radii = c(3, 30),
                                      minneigh = 2,
                                      historical = TRUE,
                                      ndays = 30,
                                      # spatial similarity arguments
                                      similarity.radii = 30,
                                      similarity.duration = 168,
                                      compl.similarity.threshold = 90,
                                      # correction arguments
                                      mdata = NULL,
                                      mdata.column = NULL,
                                      mid.column = NULL,
                                      m.radii = 30,
                                      historical.corr = TRUE,
                                      ndays.corr = 30,
                                      nmin = 5,
                                      compl.correction.threshold = 90,
                                      conf.level = 0.95,
                                      robust = TRUE,
                                      split = TRUE,
                                      # quality flag arguments
                                      qcflags = c('flag_range', 'flag_constant_value', 'flag_outlier', 'flag_correlation', 'flag_similarity'),
                                      # other arguments
                                      x.column = 'longitude',
                                      y.column = 'latitude',
                                      param = 'PM2.5',
                                      keep = c('flags', 'none')) {
  
  # reference data are required for the correction
  if (missing(rdata) || is.null(rdata)) {
    stop("'rdata' (reference observations) is required for the correction.", call. = FALSE)
  }
  stopifnot(is.data.frame(rdata))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  # argument checks
  keep <- match.arg(keep)
  #
  # temporary column with the measurements that passed the first tests
  flagged.column <- paste0(data.column, '_flagged')
  # individual quality flags of the tests
  flag.columns <- paste0('flag', c('_range', '_constant_value', '_outlier', '_correlation', '_similarity'))
  #
  # PART 1: FILTER approach up to 'Outlier Detection' (raw data)
  data <- range_validity_qctest(data,
                                data.column = data.column,
                                upper.bound = upper.bound,
                                lower.bound = lower.bound,
                                new = FALSE,
                                param = param)
  data <- constant_value_qctest(data,
                                date.column = date.column,
                                data.column = data.column,
                                id.column = id.column,
                                persist.duration = persist.duration,
                                compl.persist.duration = compl.persist.duration,
                                min.variation = min.variation,
                                method = method,
                                new = FALSE,
                                param = param)
  data <- outlier_qctest(data,
                         date.column = date.column,
                         data.column = data.column,
                         id.column = id.column,
                         outlier.duration = outlier.duration,
                         compl.outlier.duration = compl.outlier.duration,
                         outlier.threshold = outlier.threshold,
                         new = FALSE,
                         param = param)
  data <- spatial_outlier_qctest(data,
                                 date.column = date.column,
                                 data.column = data.column,
                                 out.flag.column = out.flag.column,
                                 id.column = id.column,
                                 x.column = x.column,
                                 y.column = y.column,
                                 outlier.radii = outlier.radii,
                                 outlier.minneigh = outlier.minneigh,
                                 new = FALSE,
                                 param = param)
  #
  # PART 2: Correction of the raw measurements, with the model trained on the measurements
  # that passed the range, constant value and outlier tests
  data <- data %>%
    dplyr::mutate(!!flagged.column := ifelse(.data[['flag_range']] %in% 1 & .data[['flag_constant_value']] %in% 2 & .data[['flag_outlier']] %in% 2,
                                             .data[[data.column]], NA_real_))
  data <- correction_process(data = data,
                             date.column = date.column,
                             data.column = data.column,
                             data.flag.column = flagged.column,
                             id.column = id.column,
                             rdata = rdata,
                             rdata.column = rdata.column,
                             rid.column = rid.column,
                             mdata = mdata,
                             mdata.column = mdata.column,
                             mid.column = mid.column,
                             x.column = x.column,
                             y.column = y.column,
                             m.radii = m.radii,
                             nmin = nmin,
                             compl.correction.threshold = compl.correction.threshold,
                             conf.level = conf.level,
                             robust = robust,
                             param = param,
                             ndays.corr = ndays.corr,
                             historical.corr = historical.corr,
                             split = split)
  # removal of the temporary column and of the raw-data flags (the tests are repeated on the corrected data)
  data <- data %>%
    dplyr::select(-dplyr::any_of(c(flagged.column, flag.columns)))
  #
  # PART 3: FILTER approach (Corrected Processing Level)
  # name of the corrected values given by correction_process()
  new.data.column <- paste0(param, '_cor')
  data <- runFILTER_raw_level(data,
                              date.column = date.column,
                              data.column = new.data.column,
                              id.column = id.column,
                              upper.bound = upper.bound,
                              lower.bound = lower.bound,
                              persist.duration = persist.duration,
                              compl.persist.duration = compl.persist.duration,
                              min.variation = min.variation,
                              method = method,
                              outlier.duration = outlier.duration,
                              compl.outlier.duration = compl.outlier.duration,
                              outlier.threshold = outlier.threshold,
                              out.flag.column = out.flag.column,
                              outlier.radii = outlier.radii,
                              outlier.minneigh = outlier.minneigh,
                              rdata = rdata,
                              rdata.column = rdata.column,
                              rid.column = rid.column,
                              compl.corr.threshold = compl.corr.threshold,
                              corr.radii = corr.radii,
                              minneigh = minneigh,
                              historical = historical,
                              ndays = ndays,
                              similarity.radii = similarity.radii,
                              similarity.duration = similarity.duration,
                              compl.similarity.threshold = compl.similarity.threshold,
                              qcflags = qcflags,
                              x.column = x.column,
                              y.column = y.column,
                              param = param,
                              keep = keep)
  #
  return(data)
}