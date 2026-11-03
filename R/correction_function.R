#' Correction Algorithm
#'
#' Corrects the measurements at a location using data from nearby official instruments and, optionally, auxiliary parameters (e.g., meteorological factors). Returns the corrected values (`<param>_cor`) together with their prediction intervals (`<param>_cor_lower` and `<param>_cor_upper`).
#'
#'
#' @details
#' Let \eqn{X_{it}} denote the measurements at a location \eqn{i}, \eqn{Y_{it}} the reference data from the nearby official monitoring stations, and \eqn{Z_{n,it}} optional covariates (e.g., meteorological factors). The correction is performed in four steps:
#'
#' \strong{Step 1: Initial regression model}
#'
#' An initial regression model, \eqn{Y_{it} = f(X_{it}, Z_{n,it})}, is fitted:
#' \deqn{
#' Y_{it} = \beta_0 + \beta_x X_{it} + \sum_{n} \beta_n Z_{n,it} + \varepsilon_t,
#' }{Y = b0 + bx * X + sum(bn * Zn) + e}
#' \itemize{
#'   \item \eqn{\beta_0, \beta_x, \beta_n}: regression coefficients.
#'   \item \eqn{\varepsilon_t}: random error (residual).
#' }
#'
#' The coefficients are estimated either by Ordinary Least Squares (OLS) or by robust regression (M-estimation with Huber weights), depending on `robust`:
#' \itemize{
#'   \item OLS (`robust = FALSE`):
#'   \deqn{
#'   \hat{\boldsymbol{\beta}}^{OLS} = \arg \min_{\boldsymbol{\beta}} \sum_{t} \left( Y_{it} - \hat{Y}_{it} \right)^2
#'   }{beta_OLS = argmin sum((Y - Yhat)^2)}
#'   \item Robust (`robust = TRUE`):
#'   \deqn{
#'   \hat{\boldsymbol{\beta}}^{Robust} = \arg \min_{\boldsymbol{\beta}} \sum_t \rho \!\left( \frac{Y_{it} - \hat{Y}_{it}}{\hat{\sigma}} \right)
#'   }{beta_Robust = argmin sum(rho((Y - Yhat) / sigma))}
#'   where \eqn{\rho(\cdot)} is the Huber loss function and \eqn{\hat{\sigma} = 1.4826 \cdot \operatorname{median}\left( |\varepsilon_t - \operatorname{median}(\varepsilon)| \right)} is a robust scale estimate of the residuals, used instead of their standard deviation.
#' }
#'
#' \strong{Step 2: Detection and removal of outliers}
#'
#' The residuals of the model are \eqn{\varepsilon_t = Y_{it} - \hat{Y}_{it}}. A data point is considered an outlier if:
#' \deqn{
#' \left| \varepsilon_t - \operatorname{median}(\varepsilon) \right| \geq 3 \times \text{MAD},
#' }{|e - median(e)| >= 3 * MAD}
#' where the median absolute deviation is:
#' \deqn{
#' \text{MAD} = \operatorname{median} \left( \left| \varepsilon_t - \operatorname{median}(\varepsilon) \right| \right).
#' }{MAD = median(|e - median(e)|)}
#' Only the timestamps \eqn{t^{\ast}} that are not detected as outliers are used to develop the final model.
#'
#' \strong{Step 3: Updated regression model with significant predictors}
#'
#' The model is refitted at the timestamps \eqn{t^{\ast}}, keeping only the predictors that are significant at the confidence level \eqn{1-\alpha} (`conf.level`, default `0.95`, i.e., \eqn{\alpha = 0.05}):
#' \deqn{
#' Y_{it^{\ast}} = \beta_0 + \sum_{p \in \mathcal{P}^{\ast}} \beta_p W_{p,it^{\ast}} + \varepsilon_t,
#' }{Y = b0 + sum(bp * Wp) + e}
#' where \eqn{W_p} denotes the predictors (\eqn{X_{it}} and the covariates) and \eqn{\mathcal{P}^{\ast}} is the set of significant predictors. The significance of each coefficient \eqn{\beta_p} is assessed with the t-statistic:
#' \deqn{
#' t_p = \frac{\hat{\beta}_p}{SE(\hat{\beta}_p)}, \quad df = N_{t^{\ast}} - k - 1,
#' }{t_p = beta_p / SE(beta_p), df = N - k - 1}
#' \itemize{
#'   \item \eqn{SE(\hat{\beta}_p)}: standard error of \eqn{\hat{\beta}_p}.
#'   \item \eqn{N_{t^{\ast}}}: number of timestamps \eqn{t^{\ast}} used for the model fitting.
#'   \item \eqn{k}: number of predictors in the model.
#' }
#' The corresponding \eqn{P}-values are \eqn{P_p = 2 \, P(t > |t_p|)}{P_p = 2 * P(t > |t_p|)}, and only the predictors with \eqn{P_p < \alpha} are kept in the final model.
#'
#' \strong{Step 4: Predictions with prediction intervals}
#'
#' For each timestamp \eqn{t}, the corrected values (`<param>_cor`) are calculated as:
#' \deqn{
#' \hat{Y}_{it} = \hat{\beta}_0 + \sum_{p \in \mathcal{P}^{\ast}} \hat{\beta}_p W_{p,it},
#' }{Yhat = b0 + sum(bp * Wp)}
#' and their prediction interval (`<param>_cor_lower` and `<param>_cor_upper`) as:
#' \deqn{
#' \hat{Y}_{it} \pm t_{\alpha/2,df} \cdot \sqrt{\hat{\sigma}^2 + SE(\hat{Y}_{it})^2},
#' }{Yhat +- t * sqrt(sigma^2 + SE(Yhat)^2)}
#' where \eqn{t_{\alpha/2,df}} is the critical value of the t-distribution with \eqn{df} degrees of freedom, \eqn{\hat{\sigma}^2} is the residual variance of the final model, and \eqn{SE(\hat{Y}_{it})} is the standard error of the fitted value.
#'
#' @param df A data frame containing the dataset used for model fitting and prediction.
#' @param y character. Name of the dependent variable (response) column in `df`.
#' @param x character. Name of the main independent variable (predictor) column in `df`.
#' @param covariates character. Vector of covariate column names in `df`. Can be `NULL` if no additional parameters are used.
#' @param conf.level numeric. Confidence level for the selection of significant predictors and the prediction intervals. Default is `0.95`.
#' @param param character. Base name of the output columns (`<param>_cor`, `<param>_cor_lower`, `<param>_cor_upper`). Default is `"PM2.5"`.
#' @param newdata A data frame for which predictions are also performed. Must contain the columns specified by `x` and `covariates`. Default is `NULL`.
#' @param robust logical. If `TRUE`, robust regression (`MASS::rlm`) is used; otherwise, ordinary least squares (`lm`). Default is `TRUE`.
#'
#' @return A list with two data frames:
#' \itemize{
#'   \item \code{target}: `df` with three new columns: \code{<param>_cor} (corrected values), \code{<param>_cor_lower} and \code{<param>_cor_upper} (lower and upper bounds of the prediction intervals). They are `NA` where the correction could not be calculated.
#'   \item \code{new}: if `newdata` is provided, a data frame with the same three columns for `newdata`; otherwise, `NULL`.
#' }
#'
#' @export
#'
#' @seealso \code{\link{correction_process}}
#'
#' @export
correction_function <- function(df,
                                y = 'reference',
                                x = 'sensor',
                                covariates = 'rh',
                                conf.level = 0.95,
                                newdata = NULL,
                                param = 'PM2.5',
                                robust = TRUE) {
  # argument flags
  stopifnot(is.character(y), y %in% colnames(df))
  stopifnot(is.character(x), x %in% colnames(df))
  stopifnot(is.null(covariates) ||
              (is.character(covariates) && all(covariates %in% colnames(df))))
  stopifnot(is.character(param), !is.null(param))
  stopifnot(is.numeric(conf.level), !is.null(conf.level))
  stopifnot(is.logical(robust))
  #
  # full formula
  model.formula <- if (is.null(covariates)) {
    reformulate(x, response = y)
  } else {
    reformulate(c(x, covariates), response = y)
  }
  # the warnings of rlm (e.g., 'rlm' failed to converge) are suppressed
  initial.model <- if (robust) {
    suppressWarnings(MASS::rlm(model.formula, data = df, na.action = na.exclude))
  } else {
    lm(model.formula, data = df, na.action = na.exclude)
  }
  # remove outliers using MED-MAD function
  df.filtered <- df %>%
    dplyr::mutate(Bias = residuals(initial.model)) %>%
    dplyr::filter(abs(.data[['Bias']] - median(.data[['Bias']], na.rm = TRUE)) <
                    3 * mad(.data[['Bias']], constant = 1, na.rm = TRUE)) %>%
    dplyr::select(-dplyr::all_of('Bias'))
  #
  # correction model
  if (robust) {
    model <- suppressWarnings(MASS::rlm(model.formula, data = df.filtered, na.action = na.exclude))
    model.summary <- suppressWarnings(summary(model))
    t_statistic <- coefficients(model.summary)[, 3]
    dof <- model.summary$df[2]
    p.params <- 2 * pt(abs(t_statistic), dof, lower.tail = FALSE)
  } else {
    model <- lm(model.formula, data = df.filtered, na.action = na.exclude)
    p.params <- coefficients(summary(model))[, 4]
  }
  #
  # x is always kept; covariates are kept only if significant
  sig.params <- names(p.params)[p.params < (1 - conf.level)]
  sig.covariates <- covariates[vapply(covariates, function(cv) {
    if (is.numeric(df[[cv]])) {
      cv %in% sig.params                  # numeric: coefficient name == column name
    } else {
      any(startsWith(sig.params, cv))     # factor/character: e.g. 'seasonWinter'
    }
  }, logical(1))]
  model.params <- c(x, sig.covariates)
  #
  # updated model
  model.formula <- reformulate(model.params, response = y)
  if (robust) {
    updated.model <- suppressWarnings(MASS::rlm(model.formula, data = df.filtered, na.action = na.exclude))
  } else {
    updated.model <- lm(model.formula, data = df.filtered, na.action = na.exclude)
  }
  #
  # predictions
  # the warnings of the prediction intervals are suppressed
  df.predictions <- suppressWarnings(
    predict(updated.model, newdata = df, interval = 'confidence', level = conf.level)
  ) %>% tibble::as_tibble()
  df <- df %>%
    dplyr::mutate(!!paste0(param, '_cor') := round(df.predictions[['fit']], 1),
                  !!paste0(param, '_cor_upper') := round(df.predictions[['upr']], 1),
                  !!paste0(param, '_cor_lower') := round(df.predictions[['lwr']], 1))
  if (!is.null(newdata)) {
    newdata.predictions <- suppressWarnings(
      predict(updated.model, newdata = newdata, interval = 'prediction', level = conf.level)
    ) %>% tibble::as_tibble()
    newdata <- newdata %>%
      dplyr::mutate(!!paste0(param, '_cor') := round(newdata.predictions[['fit']], 1),
                    !!paste0(param, '_cor_upper') := round(newdata.predictions[['upr']], 1),
                    !!paste0(param, '_cor_lower') := round(newdata.predictions[['lwr']], 1))
  }
  return(list(target = df, new = newdata))
}