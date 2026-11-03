#' Package imports
#'
#' @name utils
#' @import dplyr
#' @import tidyr
#' @import runner
#' @import lubridate
#' @import zoo
#' @import sf
#' @import tibble
#' @importFrom stats sd var cor reformulate lm residuals median mad coefficients pt predict na.exclude
#' @importFrom MASS rlm
#' @importFrom rlang .data := is_missing
NULL

# ------------------------------------------------------------------
# Internal helper functions (not exported)
# ------------------------------------------------------------------

Mode_fun <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) {
    if (is.factor(x)) return(factor(NA, levels = levels(x)))
    if (is.character(x)) return(NA_character_)
    if (is.logical(x)) return(NA)
    if (is.integer(x)) return(NA_integer_)
    if (is.numeric(x)) return(NA_real_)
    return(NA)
  }
  modes <- tibble(x = x) %>%
    count(x, sort = TRUE) %>%
    filter(n == max(n)) %>%
    select(x)
  if (is.factor(x)) {
    return(factor(modes, levels = levels(x)))
  } else if (is.character(x)) {
    return(as.character(modes))
  } else if (is.logical(x)) {
    return(as.logical(modes))
  } else if (is.integer(x)) {
    return(as.integer(modes))
  } else if (is.numeric(x)) {
    return(as.numeric(modes))
  } else {
    return(modes)
  }
}
