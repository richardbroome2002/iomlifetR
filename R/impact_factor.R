#' Calculate the Impact Factor in future years
#'
#' @param pm_conc A numeric vector of PM2.5 concentrations in future years. The final value
#' will be assumed to remain constant until the end of the 120 year  period.
#' @param lag_structure A numeric vector specifying the structure of any cessation lag.
#' @param RR A number. Specifies the relative risk from an epidemiologiccal study.
#' @param unit A number. The unit change associated with the relative  risk
#' @return A numeric vector of length 120 specifying the impact factor over the next 120 years
#' @description produces a vector of the impact factor over 120 years.
#' @examples
#'
#' # No cessation lag. PM concentration falls by 1 mcg/m3 after 1 year
#' pm <- c(15, 14)
#' impact_factor(pm)
#'
#' #US EPA cessation-lag structure
#' us_epa_lag <- cumsum(c(0.3, rep(0.5, 4), rep(0.2/15, 15)))
#' impact_factor(pm, us_epa_lag)
#'
#' # 5 year linear cessation-lag. PM concentration falls by 0.1mcg/m3 per year for 10 years.
#' pm <- rev(seq(15, 14, 0.1))
#' linear_lag <- seq(0.1, 1, 5)
#' impact_factor(pm, linear_lag)

impact_factor <- function(pm_conc,
                          lag_structure = 1,
                          RR = 1.062, unit = 10){

  # Create a variable that allows us to get all the variables the right length
  duration <- 120
  pm <- c(pm_conc, rep(pm_conc[length(pm_conc)], duration - length(pm_conc)))
  # Calculate the year-to-year change in PM concentration.
  change <- c(0, diff(pm))
  # Which years have a change in PM concentration?
  years <- which(change != 0)

  # Make a matrix of the PM2.5 concentration adjusted for the lag effect
  x <- matrix(0, duration, duration)
  for(i in seq_len(duration)){
    x[i, ] <- c(rep(0, i - 1),
                head(lag_structure, duration - i), # If we're getting to the end of the row the lag has to be cut.
                rep(1, length(pm) - i + 1 - length(head(lag_structure, duration - i)))
    )
  }

  x <- x * -change  # lag x the change in concentration
  x <- colSums(x)   # Total effect in each year

  # Calculate the impact factor
  IF <- 1 / RR^(x/unit)
  IF
}
