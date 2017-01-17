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
#' points(no_lag$year, no_lag$deaths_ext, col = "red")
#' abline(h = 0)
#' legend(2100, 700, c("Current cohort", "Extended population"),
#'        col=c("black", "red"), pch = 16, cex = 1.5, lty=2)
#' plot(no_lag$year, no_lag$ly_current,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", , type = "b")
#' points(no_lag$year, no_lag$ly_ext, col = "red")
#' abline(h = 0)
#'
#' # US EPA lag
#' lag <- cumsum(c(0.3, rep(0.5/4, 4), rep(0.2/15, 15)))
#' epa_lag <- impact(demog_data, lag_structure = lag)
#'
#' # Comparison of no lag and US EPA lag (Extended cohort)
#' plot(no_lag$year, no_lag$deaths_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(epa_lag$year, epa_lag$deaths_ext, col = "red", type = "b")
#' abline(h = 0)
#' legend(2100, 700, c("No lag", "US EPA lag"),
#'        col=c("black", "red"), pch = 21, cex = 2, lty=2)
#' plot(no_lag$year, no_lag$ly_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", , type = "b")
#' points(epa_lag$year, epa_lag$ly_ext, col = "red", type = "b")
#' abline(h = 0)
#'
#' # Assuming PM takes 10 years to fall by 1mcg and US EPA cessation lag
#' pm <- seq(0.1, 1, 0.1)
#' slow_pm <- impact(demog_data, delta_pm = pm, lag_structure = lag)
#'
#' plot(no_lag$year, no_lag$deaths_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Deaths avoided", type = "b")
#' points(epa_lag$year, epa_lag$deaths_ext, col = "red", type = "b")
#' points(slow_pm$year, slow_pm$deaths_ext, col = "blue", type = "b")
#' abline(h = 0)
#' legend(2080, 700, c("No lag", "US EPA lag", "Lag and gradual fall in PM"),
#'        col=c("black", "red", "blue"), pch = 21, cex = 2, lty=2)
#' plot(no_lag$year, no_lag$ly_ext,
#'      xlab = "Year", ylab = "Number",
#'      main = "Life-years gained", type = "b")
#' points(epa_lag$year, epa_lag$ly_ext, col = "red", type = "b")
#' points(slow_pm$year, slow_pm$ly_ext, col = "blue", type = "b")
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
    Mx[min_age_at_risk:length(Mx)] <-  Mx[min_age_at_risk:length(Mx)] * IF
    n <- rep(1, length(Mx))

    # ax: adjustment for timing of death within each age-group
    if(neonatal_deaths){
      ax <- c(0.1, rep(0.5, length(n) - 1))
    } else {
      ax <- rep(0.5, length(n))
    }

    qx <- n * Mx / (1 + n * (1 - ax) * Mx)

    # Calculate qx at ages older than open ended age groups assuming a log-linear increase
    mod_data <- data.frame(age = 50:(length(qx) - 2))
    mod_data$lqx <- log(qx[50:(length(qx) - 2)])
    mod <- lm(lqx ~ age, data = mod_data)
    qx[length(qx):(max_age + 1)] <- exp(predict(mod, newdata = data.frame(age = length(qx):(max_age + 1))))

    qx[length(qx)] <- 1
    # Sx: survival probability
    sx <- 1 - qx
    sx
    }

  # Make a list of Leslie matricies (ie a matrix where the diagonal is a vector
  # of age-specific survival probabilities)
  make_leslie_matrices <- function(sx_matrix) {
    # Make extra values to buffer each matrix. The top left-hand corner of each
    # matrix is 1, so that the number of births remains constant). Potentially,
    # this could be changed so that the number of births each year could be
    # varied.
    r1 <- rep(c(1,0), c(1, max_age - 1))
    c1 <- rep(0, max_age + 1)
    # Make a diagnoal matrix from  each column of the matrix of survival probabilities (minus the last row)
    sx_diag <- lapply(data.frame(sx_matrix[-nrow(sx_matrix),1:ncol(sx_matrix)]), diag)
    # Bind the buffers to each matrix
    sx_diag <- lapply(sx_diag, function(x, y){ rbind(y, x) }, y = r1)
    sx_diag <- lapply(sx_diag, cbind, c1)
    sx_diag
  }

  # Make a matrix of future population by age and calendar year
  population_matrix <- function(leslie_matrix_list) {
    # Create an empty matrix with rows representing age (0 - max age) and
    # columns representing years.
    pop_mat <- matrix(0,
                      ncol = length(leslie_matrix_list) + 1,
                      nrow = nrow(leslie_matrix_list[[1]]))
    # Make a vector of population. The open ended age group is distributed amont
    # all ages up to max_age assuming the population delines linearly
    nx <- c(demog_data$population)
    open_pop <- nx[length(nx)]
    years <- (max_age + 2) - (length(nx) - 1)
    first_year = 2 * open_pop / years
    pop <- first_year - (first_year/years) * 1:years - 1
    # pop[length(pop)] <- 0
    nx <- c(nx[-length(nx)], pop[-length(pop)])

    # This makes the first column of pop_mat.
    pop_mat[, 1] <- nx
    # We then iterate along, multiplying each column of population by the appropriate Leslie matrix
    for(i in seq_along(leslie_matrix_list)){
      pop_mat[, i + 1] <- leslie_matrix_list[[i]] %*% pop_mat[, i]
    }
    # Remove the last column and do some renaming.
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

  data.frame(
    year = base_year:(base_year + 119),
    deaths_ext = colSums(diff_deaths),
    deaths_current = colSums(diff_deaths_current),
    ly_ext = colSums(diff_ly),
    ly_current = colSums(diff_ly_current)
  )
}


# Present value -----------------------------------------------------------

#' Present value of life years gained in the future
#'
#' @param impact_obj A dataframe - the output of the `impact()` function
#' @param vly The value of a statistical life year
#' @param discount_rate The discount rate to be applied
#' @param time_horizon The number of future years to be included
#'
#' @return A list with:
#' \itemize{
#'       \item{The present value of future life years gained in the cohort alive in the first year of assessment}
#'       \item{The present value of future life years gained in a population refreshed by new births}
#'       \item{A data frame showing the present value of life yearas }
#'       }
#'
#' @export
#'
#' @examples
#'
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
#' val_3 <-  present_value(no_lag, 250000)
#' val_1 <-  present_value(no_lag, 250000, discount_rate = 0.01)
#' val_7 <-  present_value(no_lag, 250000, discount_rate = 0.07)
#'
#' plot(val_1[[3]]$year, val_1[[3]]$pv_ext,
#'      xlab = "Year", ylab = "Present value ($)",
#'      main = "Present value of life years gained")
#' points(val_3[[3]]$year, val_3[[3]]$pv_ext, col="red")
#' points(val_7[[3]]$year, val_7[[3]]$pv_ext, col="blue")
#' legend(2090, 2.7e9, c("1% discount", "3% discount", "7% discount"),
#'        col=c("black", "red", "blue"), pch = 1, cex = 1)

present_value <- function(impact_obj, vly, discount_rate = 0.03, time_horizon = 106) {

  # Add a vector of time in years to the impact data frame
  impact_obj$t <- 1:nrow(impact_obj)
  # Calculate the present value of life years saved in the current cohort
  impact_obj$pv_current <- (impact_obj$ly_current * vly) / (1 + discount_rate)^impact_obj$t
  # Calculate the present value of life years saved in the extended cohort
  impact_obj$pv_ext <- (impact_obj$ly_ext * vly) / (1 + discount_rate)^impact_obj$t

  # Prepare results
  results <- list(pv_current <- sum(impact_obj$pv_current[1:time_horizon]),
                  pv_ext <- sum(impact_obj$pv_ext[1:time_horizon]),
                  annual_results <- data.frame(year = impact_obj$year[1:time_horizon],
                                               ly_current = impact_obj$ly_current[1:time_horizon],
                                               pv_current = impact_obj$pv_current[1:time_horizon],
                                               ly_ext = impact_obj$ly_current[1:time_horizon],
                                               pv_ext = impact_obj$pv_ext[1:time_horizon]
                                               ))
  return(results)
}
