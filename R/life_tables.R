#' Calculate life expectancy
#'
#' @descripiton Calculate life expectancy using the Chiang method
#' @param start_age A numeric vector. The first age (in years) in each age group.
#' @param hazard A numeric vector. The age specific hazard.
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
#' Example using data from
#' ons_data <- read.csv("ons_data.csv")
#' start_age <- c(0, 1, seq(5, 85, 5))
#' hazard <- ons_data$deaths / ons_data$pop
#' life_expectancy(start_age, hazard)
#'
life_expectancy <- function(start_age, hazard, neonatal_deaths = TRUE){

  if (length(start_age) != length(hazard))
  {warning("The number of age groups is not equal to the number of hazards")}

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


#' An implementation of the IOM life expectancy spreadsheets
#'
#' @description Calculate the change in life expectancy associated with a reduction in risk of death
#' @param start_age Numeric vector. The starting age of each age group
#' @param population A numeric vector. The age-specific population size.
#' @param deaths A numeric vector. The age-specific number of deaths.
#' @param ages_at_risk Numeric vector. The age groups exposed to risk.
#' @param pm_concentration A number. The population weighted-mean PM2.5 concentration of interest.
#' @param RR A number. Specifies the relative risk from an epidemiologiccal study.
#' @param neonatal_deaths Logical. Are neonatal deaths included?
#' @return A numeric vector of length 120 specifying the impact factor over the next 120 years.
#'
#' @return A list of interesting outputs from the burden calculations
#' @export
#'
#' @examples
#'
#' current_burden(0:90, population, age)
#'
calculate_burden <- function(start_age, population, deaths, ages_at_risk = 30:90,
                             pm_concentration = 1, RR = 1.06, unit = 10,
                             neonatal_deaths = TRUE){
  # Caclulate additional risk associated with.
  beta   <- log(RR) / unit
  rr     <- exp(-beta * pm_concentration)

  # Align rr with ages at risk. If age is not at risk, rr = 1.
  impact <- ifelse(start_age %in% ages_at_risk, rr, 1)

  baseline_hazard <- deaths / population
  impacted_hazard <- baseline_hazard * impact

  # Produce life tables for the baseline and reduced exposure scenarios
  baseline_lt <- life_expectancy(start_age, baseline_hazard, neonatal_deaths)
  impacted_lt <- life_expectancy(start_age, impacted_hazard, neonatal_deaths)

  differences <- data.frame(
    age <- baseline_lt$age,
    ex_diff = 365 * (impacted_lt$ex - baseline_lt$ex),
    ly_diff = impacted_lt$Lx - baseline_lt$Lx,
    dx_diff = baseline_lt$dx - impacted_lt$dx
  )
  results <- list(
    baseline = baseline_lt,
    impacted = impacted_lt,
    difference = differences
  )
  return(results)
}
