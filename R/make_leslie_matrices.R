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
