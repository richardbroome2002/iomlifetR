# Impact Factor -----------------------------------------------------------

#' Calculate the Impact Factor in future years
#' @description produces a vector of the impact factor over 120 years.
#' @param delta_pm A numeric vector. The annual reduction from baseline of
#'   future PM2.5 concentrations. The final value will be assumed to remain
#'   constant until the end of the 120 year  period.
#' @param lag_structure A numeric vector specifying the structure of any
#'   cessation lag.
#' @param RR A number. Specifies the relative risk from an epidemiologiccal
#'   study.
#' @param unit A number. The unit change associated with the relative  risk
#' @return A numeric vector of length 120 specifying the impact factor over the
#'   next 120 years
#' @export
#' @examples
#' # PM concentration reduced by 1 mcg/m3 in year 1, no cessation lag, RR = 1.06 per 10mcg/m3.
#' impact_factor()
#'
#' #US EPA cessation-lag structure
#' epa_lag <- cumsum(c(0.3, rep(0.5, 4), rep(0.2/15, 15)))
#' impact_factor(lag_structure = epa_lag)
#'
#' # 5 year linear cessation-lag. PM concentration falls by 0.1mcg/m3 per year for 10 years.
#' pm <- seq(0.1, 1, 0.1)
#' linear_lag <- seq(0.1, 1, 5)
#' impact_factor(pm, linear_lag)

impact_factor <- function(delta_pm = 1, lag_structure = 1, RR = 1.06, unit = 10){

  # Create a variable that allows us to get all the variables the right length
  duration <- 120
  pm <- c(0, delta_pm, rep(delta_pm[length(delta_pm)], duration - length(delta_pm)))
  # Calculate the year-to-year change in PM concentration.
  change <- diff(pm)
  # Which years have a change in PM concentration?
  # years <- which(change != 0)

  # Make a matrix of the PM2.5 concentration adjusted for the lag effect
  x <- matrix(0, duration, duration)
  for(i in seq_len(duration)){
    x[i, ] <- c(rep(0, i - 1),
                head(lag_structure, duration - i), # If we're getting to the end of the row the lag has to be cut.
                rep(1, length(change) - i + 1 - length(head(lag_structure, duration - i)))
    )
  }
  x <- x * change  # lag times the change in concentration
  x <- colSums(x)   # Total effect in each year


  # Calculate the impact factor
  IF <- 1 / RR^(x/unit)
  IF
}

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
#'
#' # Plot the effect in the current cohort and the extended population
#' par(mfrow = c(2, 1), cex = 0.5, cex.main = 1.9,
#'     cex.lab = 1.5, cex.axis = 1.4)
#' plot(no_lag$year, no_lag$deaths_current,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(no_lag$year, no_lag$deaths, col = "red")
#' abline(h = 0)
#' legend(2100, 700, c("Current cohort", "Total"),
#'        col=c("black", "red"), pch = 16, cex = 1.5, lty=2)
#' plot(no_lag$year, no_lag$ly_current,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", , type = "b")
#' points(no_lag$year, no_lag$ly, col = "red")
#' abline(h = 0)
#'
#' # US EPA lag
#' lag <- cumsum(c(0.3, rep(0.5/4, 4), rep(0.2/15, 15)))
#' epa_lag <- impact(demog_data, lag_structure = lag)
#'
#' # Comparison of no lag and US EPA lag (Extended cohort)
#' plot(no_lag$year, no_lag$deaths,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(epa_lag$year, epa_lag$deaths, col = "red", type = "b")
#' abline(h = 0)
#' legend(2100, 700, c("No lag", "US EPA lag"),
#'        col=c("black", "red"), pch = 21, cex = 2, lty=2)
#' plot(no_lag$year, no_lag$ly,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", , type = "b")
#' points(epa_lag$year, epa_lag$ly, col = "red", type = "b")
#' abline(h = 0)
#'
#' # Assuming PM takes 10 years to fall by 1mcg and US EPA cessation lag
#' pm <- seq(0.1, 1, 0.1)
#' slow_pm <- impact(demog_data, delta_pm = pm, lag_structure = lag)
#'
#' plot(no_lag$year, no_lag$deaths,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(epa_lag$year, epa_lag$deaths, col = "red", type = "b")
#' points(slow_pm$year, slow_pm$deaths, col = "blue", type = "b")
#' abline(h = 0)
#' legend(2080, 700, c("No lag", "US EPA lag", "Lag and gradual fall in PM"),
#'        col=c("black", "red", "blue"), pch = 21, cex = 2, lty=2)
#' plot(no_lag$year, no_lag$ly,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", type = "b")
#' points(epa_lag$year, epa_lag$ly, col = "red", type = "b")
#' points(slow_pm$year, slow_pm$ly, col = "blue", type = "b")
#' abline(h = 0)
impact  <- function(demog_data,
                    delta_pm = 1,
                    lag_structure = 1,
                    RR = 1.06, unit = 10,
                    max_age = 105, base_year = 2013,
                    min_age_at_risk = 30,
                    neonatal_deaths = TRUE){

  # First estimate the impact factor
  IF <- impact_factor(delta_pm = delta_pm, lag_structure = lag_structure,
                      RR = RR, unit = unit)

  # Functions -----------------------------------------
  # Calculate the survival propbability from age 0 - max_age.
  survival_probability <- function(IF) {
    # Extend the hazard and adjust it for the IF
    Mx <- c(Mx, rep(Mx[length(Mx)], max_age + 1 - length(Mx)))
    Mx[min_age_at_risk:max_age] <-  Mx[min_age_at_risk:max_age] * IF
    n <- rep(1, length(Mx))

    # ax: adjustment for timing of death within each age-group
    if(neonatal_deaths){
      ax <- c(0.1, rep(0.5, length(n) - 1))
    } else {
      ax <- rep(0.5, length(n))
    }

    qx <- n * Mx / (1 + n * (1 - ax) * Mx)
    qx[length(qx)] <- 1
    # Sx: survival probability
    sx <- 1 - qx
    sx

    # Extend the hazard by repeating hazard in the final age group to length 106
    # Mulitply ages 30+ by the impact factor

  }

  # Make a Leslie matrix to transform the current population
  survival_diag <- function(IF){
    # Calculate the survival probability
    sx <- survival_probability(IF)
    # Make the diagonal
    sx_diag <- diag(head(sx, - 1))
    # Make a vector that goes 1, 0, 0 ... and make it the first row
    r1 <- rep(c(1,0), c(1, max_age - 1))
    sx_diag <- cbind(rbind(r1, sx_diag), rep(0, max_age + 1))
    sx_diag
  }

  # Make a matrix of future population by age and calendar year
  population_matrix <- function(IFs) {
    # Make a list of Leslie matrices
    sx_diags <- lapply(IFs, survival_diag)
    # Empty matrix to which the transformed popualtion data will be added
    pop_mat <- matrix(0, ncol = length(sx_diags) + 1, nrow = nrow(sx_diags[[1]]))
    # Baseline population (0 in age-groups older than we have data for)
    nx <- c(demog_data$population,
            rep(0, nrow(sx_diags[[1]]) - length(demog_data$population)))
    pop_mat[, 1] <- nx
    for(i in seq_along(sx_diags)){
      pop_mat[, i + 1] <- sx_diags[[i]] %*% pop_mat[, i]
    }
    pop_mat <- pop_mat[, -ncol(pop_mat)]
    colnames(pop_mat) <- base_year:(base_year + 119)
    rownames(pop_mat) <- 0:max_age
    pop_mat
  }

  # Calculate results --------------------------------------------------

  # Make matrices of the surival probability of the baseline and impacted
  # populations
  Mx <- demog_data$deaths/demog_data$population
  sx_mat <- sapply(rep(1, length(IF)), FUN = survival_probability)
  sx_mat_i <- sapply(IF, FUN = survival_probability)

  # Make matricies of the baseline and impacted populations
  baseline_pop <- population_matrix(rep(1, length(IF)))
  impacted_pop <- population_matrix(IF)

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

  data.frame(
    year = base_year:(base_year + 119),
    deaths_ext = colSums(diff_deaths),
    deaths_current = colSums(diff_deaths_current),
    ly_ext = colSums(diff_ly),
    ly_current = colSums(diff_ly_current)
  )
}
