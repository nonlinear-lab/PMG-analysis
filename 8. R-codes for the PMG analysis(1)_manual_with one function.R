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
# Step 1: to create panel data
#++++++++++++++++++++++++++++++++++++++++++++++++
# Reshape from wide to long format
panel_data <- my_data %>%
  pivot_longer(
    cols = -Time_Index,
    names_to = c(".value", "state"),
    names_pattern = "([xy])([0-9]+)"
  ) 
panel_data$state <- as.numeric(panel_data$state)

# Calculate lags and differences for the whole panel at once
pmg_data <- panel_data %>%
  group_by(state) %>%
  mutate(
    dy    = y - lag(y),
    dx    = x - lag(x),
    y_lag = lag(y),
    x_val = x        
  ) %>%
  filter(!is.na(dy) & !is.na(dx)) %>%
  ungroup()

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: pgm function
#++++++++++++++++++++++++++++++++++++++++++++++++
manual_pmg <- function(data_panel) {
  states     <- unique(data_panel$state)
  N          <- length(states)
  theta      <- 0
  tol        <- 1e-6    # Convergence tolerance
  converged  <- FALSE
  phi_list   <- numeric(N)

  for (iter in 1:1000) {
    theta_old <- theta
    num_sum <- 0
    den_sum <- 0
    
    for (i in 1:N) {
      group_df <- data_panel %>% filter(state == states[i])
      dy    <- as.matrix(group_df$dy)
      y_lag <- as.matrix(group_df$y_lag)
      x_val <- as.matrix(group_df$x_val)
      dx    <- as.matrix(group_df$dx)
      T_x <- nrow(dy)
      IT  <- diag(T_x)
      W1      <- cbind(1, dx)        # PMG
      IT      <- diag(T_x)
      H1      <- IT - W1 %*% solve(t(W1) %*% W1) %*% t(W1)
      xi      <- y_lag - x_val * theta #PGM
      phi_i <- solve(t(xi) %*% H1 %*% xi) %*% (t(xi) %*% H1 %*% dy) # PGM
      phi_i <- as.numeric(phi_i)
      phi_list[i] <- phi_i

      W2 <- cbind(1, dx)         # PGM
      H2 <- IT - W2 %*% solve(t(W2) %*% W2) %*% t(W2)
      res_i <- dy - phi_i * xi                                 # PGM only
      sigma2_i <- as.numeric(t(res_i) %*% H2 %*% res_i) / T_x # PMG only 
      term_den <- (phi_i^2 / sigma2_i) * as.numeric(t(x_val) %*% H2 %*% x_val)   # PMG
      term_num <- (phi_i / sigma2_i) * as.numeric(t(x_val) %*% H2 %*% (dy - phi_i * y_lag))  # PMG
      den_sum  <- den_sum + term_den
      num_sum  <- num_sum + term_num
    }
    theta_i <- - (1 / den_sum) * num_sum  # PMG
    theta <- theta_i
    print(theta)
    
    if (abs(theta - theta_old) < tol) {
      cat("Converged in", iter, "iterations.\n")
      converged <- TRUE
      break
    }
  }
  return(list(Shared_Theta = theta, Unique_Phis = phi_list))
}

results <- manual_pmg(pmg_data)
print(results)

average_phi <- round(mean(results$Unique_Phis), 4)
cat("Average Adjustment (Phi):", average_phi, "\n\n")