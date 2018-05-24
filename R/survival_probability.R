
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
