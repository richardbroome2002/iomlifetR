
# Difference in life expectancy ---------------------------------------------

#' An implementation of the IOM life expectancy spreadsheets
#'
#' @description Calculate the change in life expectancy associated with a reduction in risk of death
#' @param demog_data A data frame with columns of headed "age" (the age at which each age group begins), "population" (the size of the population) and "deaths" (the number of deaths in the population).
#' @param start_age A numeric vector. The starting age of each age group
#' @param min_age_at_risk Numeric vector. The lowest age susceptible to air pollution.
#' @param pm_concentration A number. The population weighted-mean PM2.5 concentration of interest.
#' @param RR A number. Specifies the relative risk from an epidemiologiccal study.
#' @param neonatal_deaths Logical. Are neonatal deaths included?
#' @param unit A number. Speficies the unit change associated with the relative risk (RR).
#'
#' @return A list of 3 elements:
#' \itemize{
#'   \item{The difference (between the baseline and impacted population) in:}
#'     \itemize{
#'       \item{Age-specific life expectancy (days)}
#'       \item{Age-specific life-years lived per 100,000}
#'       \item{Age-specific number of  deaths per 100,000}
#'       }
#'   \item{The baseline life table}
#'   \item{The impacted life table}}
#' @export
#'
#' @examples
#'
#' # Estimate the loss of life expectancy assoicated with 1mcg/m3
#' # of PM2.5 using an abridged set of population and mortality data:
#' population <- subset(abridged_data,
#'                      time == 2011 & sex == "Persons" & measure == "Population",
#'                      select = c(age, value))
#' population <- population$value
#' deaths <- subset(abridged_data,
#'                  time == 2011 & sex == "Persons" & measure == "Deaths",
#'                  select = c(age, value))
#' start_age <- as.numeric(gsub(" .+", "", deaths$age))
#' deaths <- deaths$value
#' demog_data <- data.frame(age = start_age, population, deaths)
#' x <- burden_le(demog_data)
#' x[[1]][1, 2] # Change in LE at birth (days)
burden_le <- function(demog_data, min_age_at_risk = 30,
                      pm_concentration = 1, RR = 1.06, unit = 10,
                      neonatal_deaths = TRUE){
  # Caclulate additional risk associated with.
  rr <- RR^(-pm_concentration / unit)
  ages_at_risk <- min_age_at_risk:max(demog_data$age)


  # Align rr with ages at risk. If age is not at risk, rr = 1.
  impact <- ifelse(demog_data$age %in% ages_at_risk, rr, 1)

  baseline_hazard <- demog_data$deaths / demog_data$population
  impacted_hazard <- baseline_hazard * impact

  # Produce life tables for the baseline and reduced exposure scenarios
  baseline_lt <- life_table(baseline_hazard, demog_data$age, neonatal_deaths)
  impacted_lt <- life_table(impacted_hazard, demog_data$age, neonatal_deaths)

  differences <- data.frame(
    age = baseline_lt$age,
    ex_diff = 365 * (impacted_lt$ex - baseline_lt$ex),
    ly_diff = impacted_lt$Lx - baseline_lt$Lx,
    dx_diff = baseline_lt$dx - impacted_lt$dx
  )
  results <- list(
    difference = differences,
    baseline = baseline_lt,
    impacted = impacted_lt
  )
  return(results)
}
