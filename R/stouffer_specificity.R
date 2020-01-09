stouffer_specificity <- function(expression_data) {
  # stouffer method:
  # (((x - x_hat) / sd(x)) + ((y - y_hat) / sd(y))) / sqrt(2)
  # where x is row and y is column
  
  means_rows <- rowMeans(expression_data)
  means_cols <- colMeans(expression_data)
  
  sd_rows <- apply(expression_data, 1, sd)
  sd_cols <- apply(expression_data, 2, sd)
  
  z_rows <- apply(expression_data, 2, function(t) (t - means_rows) / sd_rows)
  z_cols <- apply(expression_data, 1, function(t) (t - means_cols) / sd_cols)
  z_cols <- t(z_cols)
  
  z_tot <- (z_rows + z_cols) / sqrt(2)
  return(z_tot) 
}

