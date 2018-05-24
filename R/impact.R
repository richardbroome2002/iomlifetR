# Calculate Impact --------------------------------------------------------

#' An implementation of the IOMLIFET impact spreadsheets
#' @description Implements a life table method for estimating the mortality
#'   benefits of reduction in exposure to risk
#' @param demog_data A data frame with three numeric columns headed "age" (the age at which each age group begins),
#'   "population" (the size of the population) and "deaths" (the number of deaths in the population).
#' @param delta_pm The reduction in population-weighted PM2.5 concentration in
#'   each future year.
#' @param lag_structure A numeric vector with values between 0 and 1 that defines the structure of any cessation lag.
#' @param RR The relative risk (or hazard ratio) to use from the assessment.
#'   This is taken from epidemiological studies
#' @param unit The unit change in PM2.5 concentration for the RR from an epidemiological study
#' @param max_age The maximum age to use for the assessment
#' @param base_year The base year for the assessment
#' @param neonatal_deaths Logical. Are neonatal deaths included?
#'
#' @return A dataframe of five columns of:
#'     \itemize{
#'       \item{The year}
#'       \item{The difference in number of deaths over 120 years (extended population)}
#'       \item{The difference in number of deaths over 120 years (current cohort)}
#'       \item{The difference in number of Life-years over 120 years (extended population)}
#'       \item{The difference in number of Life-years over 120 years (current cohort).}
#'       }
#' @export
#'
#' @examples
#'
#' # Reduce PM by 1mcg/m3. No cessation lag. RR = 1.06 per 10mcg/m3
#' head(single_year_data)
#' population <- subset(single_year_data,
#'                      time == 2011 & sex == "Persons" & measure == "Population")
#' population <- population[, c("age", "value")]
#' population$age <- as.numeric(gsub(" .+", "", population$age))
#' colnames(population)[2] <- "population"
#' deaths <- subset(single_year_data,
#'                  time == 2011 & sex == "Persons"& measure == "Deaths")
#' deaths <- deaths[, "value"]
#' demog_data <- data.frame(population, deaths = deaths)
#'
#' no_lag <- impact(demog_data)
#' no_lag <- data.frame(lapply(no_lag, colSums))
#'
#' # Plot the effect in the current cohort and the extended population
#' par(mfrow = c(2, 1), cex = 0.5, cex.main = 1.9,
#'     cex.lab = 1.5, cex.axis = 1.4)
#' plot(rownames(no_lag), no_lag$deaths_current,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(rownames(no_lag), no_lag$deaths_ext, col = "red")
#' abline(h = 0)
#' legend(2100, 700, c("Current cohort", "Extended population"),
#'        col=c("black", "red"), pch = 16, cex = 1.5, lty=2)
#' plot(rownames(no_lag), no_lag$ly_current,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", , type = "b")
#' points(rownames(no_lag), no_lag$ly_ext, col = "red")
#' abline(h = 0)
#'
#' # US EPA lag
#' lag <- cumsum(c(0.3, rep(0.5/4, 4), rep(0.2/15, 15)))
#' epa_lag <- impact(demog_data, lag_structure = lag)
#' epa_lag <- data.frame(lapply(epa_lag, colSums))
#'
#' # Comparison of no lag and US EPA lag (Extended cohort)
#' plot(rownames(no_lag), no_lag$deaths_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(rownames(epa_lag), epa_lag$deaths_ext, col = "red", type = "b")
#' abline(h = 0)
#' legend(2100, 700, c("No lag", "US EPA lag"),
#'        col=c("black", "red"), pch = 21, cex = 2, lty=2)
#' plot(rownames(no_lag), no_lag$ly_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", , type = "b")
#' points(rownames(epa_lag), epa_lag$ly_ext, col = "red", type = "b")
#' abline(h = 0)
#'
#' # Assuming PM takes 10 years to fall by 1mcg and US EPA cessation lag
#' pm <- seq(0.1, 1, 0.1)
#' slow_pm <- impact(demog_data, delta_pm = pm, lag_structure = lag)
#' slow_pm <- data.frame(lapply(slow_pm, colSums))
#'
#' plot(rownames(no_lag), no_lag$deaths_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(rownames(epa_lag), epa_lag$deaths_ext, col = "red", type = "b")
#' points(rownames(slow_pm), slow_pm$deaths_ext, col = "blue", type = "b")
#' abline(h = 0)
#' legend(2080, 700, c("No lag", "US EPA lag", "Lag and gradual fall in PM"),
#'        col=c("black", "red", "blue"), pch = 21, cex = 2, lty=2)
#' plot(rownames(no_lag), no_lag$ly_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", type = "b")
#' points(rownames(epa_lag), epa_lag$ly_ext, col = "red", type = "b")
#' points(rownames(slow_pm), slow_pm$ly_ext, col = "blue", type = "b")
#' abline(h = 0)

impact <- function(demog_data,
                   delta_pm = 1,
                   lag_structure = 1,
                   RR = 1.06, unit = 10,
                   max_age = 105, base_year = 2013,
                   min_age_at_risk = 30,
                   neonatal_deaths = TRUE){

  # First estimate the impact factor
  IF <- impact_factor(delta_pm = delta_pm, lag_structure = lag_structure,
                      RR = RR, unit = unit)

  # Calculate results --------------------------------------------------

  # Make matrices of the surival probability of the baseline and impacted
  # populations
  Mx <- demog_data$deaths/demog_data$population
  sx_mat <- sapply(rep(1, length(IF)), FUN = survival_probability)
  sx_mat_i <- sapply(IF, FUN = survival_probability)

  # Make the leslie matricies to form the future population under the two scenarios
  baseline_leslie_matrices <- make_leslie_matrices(sx_mat)
  impacted_leslie_matricies <- make_leslie_matrices(sx_mat_i)

  # Make matricies of the baseline and impacted populations
  baseline_pop <- population_matrix(baseline_leslie_matrices)
  impacted_pop <- population_matrix(impacted_leslie_matricies)

  # Calculate number of deaths
  baseline_deaths <- baseline_pop * (1 - sx_mat)
  impacted_deaths <- impacted_pop * (1 - sx_mat_i)
  diff_deaths <- baseline_deaths - impacted_deaths

  # Calculate number of deaths among current cohort
  diff_deaths_current <- diff_deaths
  diff_deaths_current[upper.tri(diff_deaths_current)] <- 0

  # Calculate number of Life-years (population x half the number of deaths in a
  # given year)
  baseline_ly <- baseline_pop - 0.5 * baseline_deaths
  impacted_ly <- impacted_pop - 0.5 * impacted_deaths
  diff_ly <- impacted_ly - baseline_ly

  # And Life-years among the current cohort
  diff_ly_current <- diff_ly
  diff_ly_current[upper.tri(diff_ly_current)] <- 0
  results <- list(
    ly_extended = diff_ly,
    deaths_extended = diff_deaths,
    ly_current = diff_ly_current,
    deaths_current = diff_deaths_current
  )
}
