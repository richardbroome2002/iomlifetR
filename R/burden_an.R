
#' The attributable number
#'
#' @param demog_data A data frame with columns of headed "age" (the age at which each age group begins), "population" (the size of the population) and "deaths" (the number of deaths in the population).
#' @param start_age A numeric vector. The starting age of each age group
#' @param ages_at_risk A numeric vector. The age groups exposed to risk.
#' @param pm_concentration A number. The population weighted-mean PM2.5 concentration of interest.
#' @param RR A number. Specifies the relative risk from an epidemiologiccal study.
#' @param unit A number. Speficies the unit change associated with the relative risk (RR).
#'
#' @return A numeric vector of age-specific attributable numbers.
#' @export
#'
#' @examples
#' # Estimate the number of premature deaths attributable to 1mcg/m3 of PM2.5
#' data(singe_year_data)
#' year <- 2011
#' s <- "Persons"
#' population <- subset(abridged_data,
#'                      time == year & sex == s & measure == "Population",
#'                      select = c(age, value))
#' deaths <- subset(abridged_data,
#'                  time == year & sex == s & measure == "Deaths",
#'                  select = c(age, value))
#' start_age <- as.numeric(gsub(" .+", "", deaths$age))
#' demog_data <- data.frame(age = start_age, population, deaths)
#'
#' burden_an(demog_data)
burden_an <- function(demog_data, min_age_at_risk = 30,
                      pm_concentration = 1, RR = 1.06, unit = 10){
  # Calculate the RR associated with pm_concentration
  rr <- RR^(-pm_concentration / unit)
  ages_at_risk <- min_age_at_risk:max(demog_data$age)

  # Align rr with ages at risk. If age is not at risk, rr = 1.
  impact <- ifelse(demog_data$age %in% ages_at_risk, rr, 1)

  hazard <- demog_data$deaths / demog_data$population
  an <- hazard * (1 - impact) * demog_data$population
  an
}
