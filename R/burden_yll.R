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
