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
