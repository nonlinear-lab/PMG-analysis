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
# Step 1: pgm function
#++++++++++++++++++++++++++++++++++++++++++++++++
manual_pgm <- function(y, x) {
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
  
  theta   <- 0
  tol <- 1e-6          # Convergence tolerance
  converged <- FALSE
  
  for (iter in 1:1000) {
    theta_old <- theta
    num_sum <- 0
    den_sum <- 0
    T_x     <- nrow(dy)
    W1      <- cbind(1, dx)        # PMG
    IT      <- diag(T_x)
    H1      <- IT - W1 %*% solve(t(W1) %*% W1) %*% t(W1)
    xi      <- y_lag - x * theta #PGM

    phi_i <- solve(t(xi) %*% H1 %*% xi) %*% (t(xi) %*% H1 %*% dy) # PGM

    W2 <- cbind(1, dx)         # PGM
    H2 <- IT - W2 %*% solve(t(W2) %*% W2) %*% t(W2)
    res_i <- dy - phi_i * xi                                 # PGM only
    sigma2_i <- as.numeric(t(res_i) %*% H2 %*% res_i) / T_x # PMG only 
    term_den <- (phi_i^2 / sigma2_i) * as.numeric(t(x_lag) %*% H2 %*% x_lag)   # PMG
    term_num <- (phi_i / sigma2_i) * as.numeric(t(x_lag) %*% H2 %*% (dy - phi_i * y_lag))  # PMG
    den_sum  <- den_sum + term_den
    num_sum  <- num_sum + term_num
    theta_i <- - (1 / den_sum) * num_sum  # PMG
    
    theta <- theta_i
    print(theta)
    
    if (abs(theta - theta_old) < tol) {
      cat("Converged in", iter, "iterations.\n")
      converged <- TRUE
      break
    }
  }
  return(c(phi_i, theta_i))
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: run the function
#++++++++++++++++++++++++++++++++++++++++++++++++

for (j in 1:5){
  x         <- my_data[(j-1)*2+2]
  y         <- my_data[(j-1)*2+3]
  result    <- manual_pgm(y,x)
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