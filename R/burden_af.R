# Attributable Fraction ---------------------------------------------------

#' Attributable fraction
#'
#' @param RR A number. The relative risk from an epidemiological study
#' @param unit A number. The unit change in PM2.5 concentration related to RR
#' @param pm_concentration A number. The change in PM2.5 of interest.
#'
#' @return A number. The fraction of total mortality attributable to the risk factor
#' @export
#'
#' @examples

#' burden_af(RR = 1.14, pm_concentration = 0.5)
burden_af <- function(RR = 1.062, unit = 10, pm_concentration = 1) {
  rr <- RR^(pm_concentration / unit)
  af <- (rr-1)/rr
  af
}
