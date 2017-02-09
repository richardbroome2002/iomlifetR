
# Life table ---------------------------------------------------------

#' Calculate life expectancy
#'
#' @description Calculate life expectancy using the Chiang method
#' @param hazard A numeric vector. The age specific hazard.
#' @param start_age A numeric vector. The first age (in years) in each age group.
#' @param neonatal_deaths Logical. Are neonatal deaths included?
#' @return A data frame of: \itemize{
#' \item{Age groups}
#' \item{Age-specific hazards}
#' \item{Population surviving from a cohort of 100,000}
#' \item{Survival probability}
#' \item{Number of deaths per 100,000}
#' \item{Number of life-years lived per 100,0000}
#' \item{life expectancy}}
#' @export
#'
#' @examples
#'
#' # Using an abridged set of population and mortality data:
#' head(abridged_data)
#' population <- subset(abridged_data,
#'                      time == 2011 & sex == "Persons" & measure == "Population",
#'                      select = c(age, value))
#' population = population$value
#' deaths <- subset(abridged_data,
#'                  time == 2011 & sex == "Persons" & measure == "Deaths",
#'                  select = c(age, value))
#' start_age <- as.numeric(gsub(" .+", "", deaths$age))
#' deaths  <- deaths$value
#' hazard <- deaths / population
#' life_table(hazard, start_age)
#'
#' # Using single-year population and mortality data
#' head(single_year_data)
#' population <- subset(single_year_data,
#'                      time == 2011 & sex == "Persons" & measure == "Population",
#'                      select = c(age, value))
#' population = population$value
#' deaths <- subset(single_year_data,
#'                  time == 2011 & sex == "Persons" & measure == "Deaths",
#'                  select = c(age, value))
#' start_age <- as.numeric(gsub(" .+", "", deaths$age))
#' deaths  <- deaths$value
#' hazard <- deaths / population
#' life_table(hazard, start_age)
life_table <- function(hazard, start_age, neonatal_deaths = TRUE){

  # n: time within each age group
  n <- c(diff(start_age), 0)

  # ax: adjustment for timing of death within each age-group
  if(neonatal_deaths){
    ax <- c(0.1, rep(0.5, length(n) - 1))
  } else {
    ax <- rep(0.5, length(n))
  }

  # qx: mortality rate
  qx <- n * hazard / (1 + n * (1 - ax) * hazard)
  qx[length(qx)] <- 1
  # Sx: survival probability
  Sx <- 1 - qx
  # Ix: Number remaining in a cohort of 100,000
  Ix <- c(100000, 100000 * cumprod(Sx[-length(Sx)]))
  Ix_lag <- c(Ix[2:length(Ix)], 0)
  # dx: number of deaths in each age group of the cohort
  dx <- Ix - Ix_lag
  # Lx: life years lived by member so the age group
  Lx <- n * (Ix_lag + (ax * dx))
  Lx[length(Lx)] <- Ix[length(Ix)] / hazard[length(hazard)]
  # Tx: Cumulative number of life years to be lived
  Tx <- rev(cumsum(rev(Lx)))
  # ex: Life expectancy
  ex <- Tx/ Ix

  age <- ifelse(n <= 1,
                start_age,
                paste(start_age, "-", start_age + n - 1))
  age[length(age)] <- paste0(age[length(age)], "+")

  lt <- data.frame(age, hazard, Ix, Sx, dx, Lx, ex)
  lt
}


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

# Attributable number -----------------------------------------------------

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

# YLL ---------------------------------------------------------------------

#' Calculate the number of YLL
#'
#' @param attributable_number A numeric vector of age-specific attributable numbers
#' @param life_expectancy A numeriv vector of age-specific life expectancies
#' @return A numeric vector of age_specific YLL
#' @export
#'
#' @examples
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
#' le <- burden_le(demog_data)
#' impacted_le <- le[["impacted"]][, "ex"]
#' an <- burden_an(demog_data)
#' yll <- burden_yll(an, impacted_le)
#' sum(yll)
burden_yll <- function(attributable_number, life_expectancy){
  yll <- attributable_number * life_expectancy
  yll
}

# Attributable Fraction ---------------------------------------------------

#' Attributbale fraction
#'
#' @description Calculate the attributable fraction related to a PM$_{2.5}$ concentration
#' @param RR A number. Specifies the relative risk from an epidemiologiccal study.
#' @param unit A number. Speficies the unit change associated with the relative risk (RR).
#' @param pm_concentration A number. The population weighted-mean PM2.5 concentration of interest.
#'
#' @return A number. The percent of mortality attributable to the risk factor
#' @export
#'
#' @examples
#'
#' burden_af(RR = 1.14, pm_concentration = 0.5)
burden_af <- function(RR = 1.062, unit = 10, pm_concentration = 1) {
  rr <- RR^(pm_concentration / unit)
  af <- 100*(rr-1)/rr
  af
}
