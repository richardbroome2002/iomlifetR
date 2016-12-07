

# Function to calculate impact
#' Calculate the impact of a reduction in PM concentration
#'
#' @param population A vector of the population in single-year age groups
#' @param deaths A vector of the total number of deaths in the population
#' @param pm_concentration The population-weighted PM2.5 concentration in each future year
#' @param lag_structure The structure of any cessation lag
#' @param RR The relative risk (or hazard ratio) to use from the assessment. This is taken
#' from epidemiological studies
#' @param unit The unit change in PM2.5 concentration for the RR
#' @param max_age The maximum age to use for the assessment
#' @param base_year The base year for the assessment
#'
#' @return A list of four matrices (with dimensions of age and calendar year) of:
#' \itemize{
#' \item{The difference in number of deaths over 120 years (extended population)}
#' \item{The difference in number of deaths over 120 years (current cohort)}
#' \item{The difference in number of life years over 120 years (extended population)}
#' \item{The difference in number of life years over 120 years (current cohort)}}
#' @export
#'
#' @examples
#'
#' Reduce PM from 15 to 14 mcg/m3
#' pm <- c(15, 14)
#' epa_lag <- cumsum(c(0.3, rep(0.125, 4), rep(0.2/15, 15)))
#'
#' impact(pop, deaths)
impact  <- function(population, deaths,
                    pm_concentration = c(20, 19),
                    lag_structure = 1,
                    RR = 1.06, unit = 10,
                    max_age = 105, base_year = 2013){

  # First estimate the impact factor
  IF <- impact_factor(pm_conc = pm_concentration, lag_structure = lag_structure,
                      RR = RR, unit = unit)

  # Functions -----------------------------------------
  # Calculate the survival propbability from age 0 - 105.
  survival_probability <- function(IF) {
    # Extend the hazard by repeating hazard in the final age group to length 106
    qx <- c(qx, rep(qx[length(qx)], max_age + 1 - length(qx)))
    # Mulitply ages 30+ by the impact factor
    qx[31:length(qx)] <- qx[31:length(qx)] * IF
    # Calculate the surival probability in the IOMLIFET way
    sx <- (2 - qx) / (2 + qx)
    sx[length(sx)] <- 0
    sx
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
    # Baseline populaiton (0 in age-groups older than we have data for)
    nx <- c(population, rep(0, nrow(sx_diags[[1]]) - length(population)))
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
  qx <- deaths/population
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

  # Calculate number of life years (population x half the number of deaths in a
  # given year)
  baseline_ly <- baseline_pop - 0.5 * baseline_deaths
  impacted_ly <- impacted_pop - 0.5 * impacted_deaths
  diff_ly <- impacted_ly - baseline_ly

  # And life years among the current cohort
  diff_ly_current <- diff_ly
  diff_ly_current[upper.tri(diff_ly_current)] <- 0

  results <- list("Difference in deaths" = diff_deaths,
                  "Difference in deaths among current cohort" = diff_deaths_current,
                  "Difference in life years" = diff_ly,
                  "Difference in life years among current cohort" = diff_ly_current)
}
