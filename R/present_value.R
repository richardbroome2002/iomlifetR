
# Present value -----------------------------------------------------------

#' Present value of life years gained in the future
#'
#' @param life_years A numeric vector of life years gained by year
#' @param vly The value of a statistical life year
#' @param discount_rate The discount rate to be applied
#'
#' @return A numeric vector of the present value of life years gained in future years.
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
#' ly_extended <- colSums(no_lag$ly_extended)
#'
#' sum(present_value(ly_extended, 250000))
#'
present_value <- function(life_years, vly, discount_rate = 0.03) {

  # Add a vector of time in years to the impact data frame
  t <- 1:length(life_years)
  # Calculate the present value of life years saved in the current cohort
  pv <- (life_years * vly) / (1 + discount_rate)^t
  return(pv)
}
