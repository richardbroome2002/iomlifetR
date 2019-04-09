# Make a matrix of future population by age and calendar year
population_matrix <- function(leslie_matrix_list, demog_data, max_age, base_year) {
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
