#+++++++++++++++++++++++++++++++++++++++
# Step 0: Generate Sample Panel Data
#+++++++++++++++++++++++++++++++++++++++
setwd("C:/R/ardl")

# Install the ardlverse package if you haven't already
# install.packages("ardlverse")

library(tidyr)      # For reshaping the data
library(dplyr)      # For data manipulation

N <- 5    # Number of groups
T <- 200  # Number of time periods

my_data <- read.csv("data(5).csv")
mout    <- matrix(NA, nrow=5, ncol=2)

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: gm function
#++++++++++++++++++++++++++++++++++++++++++++++++
manual_gm <- function(y, x) {
  y               <- as.matrix(y)
  N_total         <- nrow(y)
  x               <- as.matrix(x)
  dx              <- diff(x)
  dy              <- diff(y)
  y_lag           <- y[-N_total]
  x_lag           <- x[-N_total]
  t_index         <- 2:N_total
  y               <- y[c(2:N_total)]
  x               <- x[c(2:N_total)]

  data <- data.frame(dy = dy,y_lag = y_lag,x_lag = x_lag,dx = dx, t_index = t_index)
  
  T_x   <- nrow(dy)
  W1   <- cbind(1, x_lag, dx)
  IT   <- diag(T_x)
  H1   <- IT - W1 %*% solve(t(W1) %*% W1) %*% t(W1)
  # xi <- y_lag - X * theta PGM
  xi <- y_lag
  
  # step 1: estimation of phi (adjustment term)
# phi_i <- solve(t(xi) %*% H_i %*% xi) %*% (t(xi) %*% H_i %*% dy) PGM
  phi_i   <- solve(t(xi) %*% H1 %*% xi) %*% (t(xi) %*% H1 %*% dy) # MG
  
  # step 2: estimation of theta (long_run term)
  W2 <- cbind(1, y_lag, dx)
  H2 <- IT - W2 %*% solve(t(W2) %*% W2) %*% t(W2)
  
  # sigma2_i <- as.numeric(t(res_i) %*% H_i %*% res_i) / T_i PMG
  # Denominator term: (phi_i^2 / sigma^2_i) * X_i' H_i X_i
  # term_den <- (phi_i^2 / sigma2_i) * as.numeric(t(X) %*% H_i %*% X)   PMG
  term_den <- phi_i * as.numeric(t(x_lag) %*% H2 %*% x_lag)           # MG
  # Numerator term: (phi_i / sigma^2_i) * X_i' H_i (dy_i - phi_i * y_{i,t-1})
  # term_num <- (phi_i / sigma2_i) * as.numeric(t(X) %*% H_i %*% (dy - phi_i * y_lag))  PMG
  term_num <- as.numeric(t(x_lag) %*% H2 %*% (dy - phi_i * y_lag))                      # MG   
  
  # theta   <- - (1 / den_sum) * num_sum    PMG
  theta_i <- - (1 / term_den) * term_num  # MG
  
  return(c(phi_i, theta_i))
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: run the function
#++++++++++++++++++++++++++++++++++++++++++++++++

for (j in 1:5){
  x         <- my_data[(j-1)*2+2]
  y         <- my_data[(j-1)*2+3]
  result    <- manual_gm(y,x)
  phi       <- result[1]
  mout[j,1] <- phi
  theta     <- result[2]
  mout[j,2] <- theta
}
print(mout)
average_theta <- round(mean(mout[, 2], na.rm = TRUE), 4)
average_phi   <- round(mean(mout[, 1], na.rm = TRUE), 4)

cat("Average Long-Run (Theta):", average_theta, "\n")
cat("Average Adjustment (Phi):", average_phi, "\n\n")