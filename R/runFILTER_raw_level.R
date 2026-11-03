#' Quality Control based on FILTER (Raw Processing Level)
#'
#' Applies the full quality control (QC) procedure of FILTER on the raw data.
#'
#' @details
#' This function performs all quality control (QC) algorithms included in the FILTER approach. QC is applied using a five-step procedure in the following order: \cr\cr
#'
#' \strong{Raw Data} ---------  \strong{Range Validity} --> \strong{Constant Value} --> \strong{Outlier Detection} --> \strong{Spatial Correlation} --> \strong{Spatial Similarity} ---------  \strong{QC Raw Data}\cr\cr
#'
#' The range validity, constant value and outlier detection (temporal and spatial) tests are applied directly on the raw data. The spatial correlation test is applied
#' only on the measurements that passed these tests (\code{flag_range = 1}, \code{flag_constant_value = 2} and \code{flag_outlier = 2}), and the spatial similarity test
#' only on the measurements that additionally passed the spatial correlation test (\code{flag_correlation = 3}). The measurements that did not pass a previous step are
#' set to \code{NA} before the next step, so that they are neither tested nor used as neighbors. At the end of the process, a total quality flag is assigned to each
#' measurement in the dataset.
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
#' @param rdata see in \code{\link{spatial_correlation_qctest}} and \code{\link{spatial_similarity_qctest}}.
#' @param rdata.column see in \code{\link{spatial_correlation_qctest}} and \code{\link{spatial_similarity_qctest}}.
#' @param rid.column see in \code{\link{spatial_correlation_qctest}} and \code{\link{spatial_similarity_qctest}}.
#' @param compl.corr.threshold see in \code{\link{spatial_correlation_qctest}}.
#' @param corr.radii see in \code{\link{spatial_correlation_qctest}}.
#' @param minneigh see in \code{\link{spatial_correlation_qctest}}.
#' @param historical see in \code{\link{spatial_correlation_qctest}}.
#' @param ndays see in \code{\link{spatial_correlation_qctest}}.
#' @param similarity.radii see in \code{\link{spatial_similarity_qctest}}.
#' @param similarity.duration see in \code{\link{spatial_similarity_qctest}}.
#' @param compl.similarity.threshold see in \code{\link{spatial_similarity_qctest}}.
#' @param qcflags see in \code{\link{quality_flag_function}}.
#' @param x.column see in \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, and \code{\link{spatial_similarity_qctest}}.
#' @param y.column see in \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, and \code{\link{spatial_similarity_qctest}}.
#' @param param see in \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, and \code{\link{spatial_similarity_qctest}}.
#' @param keep character. Which quality flags are kept in the output. If `"flags"` (default), the individual quality flags of each test are kept together with the total quality flag (`quality_flag`). If `"none"`, the individual quality flags are removed, so only the total quality flag is added.
#'
#' @return Returns the original data with the total quality flag (`quality_flag`) added as a new column for each measurement. The flag has three levels: `0` (other quality), `1` (good quality), `2` (high quality). If `keep = "flags"`, the individual quality flags of each test are also included.
#'
#' @seealso \code{\link{range_validity_qctest}}, \code{\link{constant_value_qctest}}, \code{\link{outlier_qctest}}, \code{\link{spatial_outlier_qctest}}, \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{quality_flag_function}}
#'
#' @export
runFILTER_raw_level <- function(data,
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
                                rdata = NULL,
                                rdata.column = NULL,
                                rid.column = NULL,
                                compl.corr.threshold = 90,
                                corr.radii = c(3, 30),
                                minneigh = 2,
                                historical = TRUE,
                                ndays = 30,
                                # spatial similarity arguments
                                similarity.radii = 30,
                                similarity.duration = 168,
                                compl.similarity.threshold = 90,
                                # quality flag arguments
                                qcflags = c('flag_range', 'flag_constant_value', 'flag_outlier', 'flag_correlation', 'flag_similarity'),
                                # other arguments
                                x.column = 'longitude',
                                y.column = 'latitude',
                                param = "PM2.5",
                                keep = c('flags', 'none')) {
  # argument checks
  keep <- match.arg(keep)
  #
  # temporary columns with the measurements that passed the previous steps
  flagged.column <- paste0(data.column, '_flagged')
  correlation.column <- paste0(data.column, '_correlation')
  #
  # FILTER approach (Raw Processing Level)
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
  # measurements valid in the range, constant value and outlier tests
  data <- data %>%
    mutate(!!flagged.column := ifelse(.data[['flag_range']] == 1 & .data[['flag_constant_value']] == 2 & .data[['flag_outlier']] == 2,
                                      .data[[data.column]], NA))
  data <- spatial_correlation_qctest(data,
                                     date.column = date.column,
                                     data.column = flagged.column,
                                     id.column = id.column,
                                     rdata = rdata,
                                     rdata.column = rdata.column,
                                     rid.column = rid.column,
                                     x.column = x.column,
                                     y.column = y.column,
                                     compl.corr.threshold = compl.corr.threshold,
                                     corr.radii = corr.radii,
                                     minneigh = minneigh,
                                     historical = historical,
                                     ndays = ndays,
                                     new = FALSE,
                                     param = param)
  # measurements valid also in the spatial correlation test
  data <- data %>%
    mutate(!!correlation.column := ifelse(!is.na(.data[[flagged.column]]) & .data[['flag_correlation']] == 3,
                                          .data[[data.column]], NA))
  data <- spatial_similarity_qctest(data,
                                    date.column = date.column,
                                    data.column = correlation.column,
                                    id.column = id.column,
                                    rdata = rdata,
                                    rdata.column = rdata.column,
                                    rid.column = rid.column,
                                    x.column = x.column,
                                    y.column = y.column,
                                    similarity.radii = similarity.radii,
                                    similarity.duration = similarity.duration,
                                    compl.similarity.threshold = compl.similarity.threshold,
                                    new = FALSE,
                                    param = param)
  # removal of the temporary columns
  data <- data %>%
    select(-all_of(c(flagged.column, correlation.column)))
  #
  # total quality flag
  data <- quality_flag_function(data, data.column = data.column, qcflags = qcflags)
  #
  # individual quality flags of the tests
  flag.columns <- paste0('flag', c('_range', '_constant_value', '_outlier', '_correlation', '_similarity'))
  if (keep == 'none') {
    # keep only the total quality flag
    data <- data %>%
      select(-any_of(flag.columns))
  }
  #
  return(data)
}
