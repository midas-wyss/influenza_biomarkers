# biogps data ----
load(file = "~/work/data/2019-02-13_biogps.RData")

# funtion to use ----
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

# average of tissues ----
utis <- unique(biogps_data$tissues)
avg_E <- matrix(0, nrow = nrow(biogps_data$expression), ncol = length(utis))
for (i in seq(1, length(utis))) {
  u1 <- which(biogps_data$tissues == utis[i])
  if (length(u1) > 1) {
    avg_E[, i] <- rowMeans(biogps_data$expression[, u1])
  } else {
    avg_E[, i] <- biogps_data$expression[, u1]
  }
}

# lung specific genes ----
stouf_thr <- 3
ts_stouffer <- stouffer_specificity(avg_E) # <-- see how i calculate the average of tissue expression
spec_mat <- matrix(0, nrow = nrow(x = ts_stouffer), ncol = ncol(ts_stouffer))
spec_mat[ts_stouffer >= stouf_thr] <- 1

spec_genes <- which(rowSums(spec_mat) <= 2)
lung_genes <- biogps_data$genes$SYMBOL[intersect(which(spec_mat[, 2] == 1), spec_genes)]


a1 <- spec_mat[, grep("lung", biogps_data$tissues)]
a2 <- which(rowSums(spec_mat) == 1)
intersect(a1, a2)
