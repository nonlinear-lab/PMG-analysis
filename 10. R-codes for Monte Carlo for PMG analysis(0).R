#+++++++++++++++++++++++++++++++++++
# Step 1: set monte carlo simulation
#+++++++++++++++++++++++++++++++++++
# Monte Carlo demonstration of SUR efficiency gain
library(tidyr)
library(dplyr)
set.seed(123)
reps      <- 1000
T_obs     <- 200
rho       <- 0.8
true_beta <- 2.5
T_obs     <- T
mout1 <- matrix(NA, nrow = reps, ncol = 5)
mout2 <- matrix(NA, nrow = reps, ncol = 5)
mout3 <- matrix(NA, nrow = reps, ncol = 1)
mout4 <- matrix(NA, nrow = reps, ncol = 1)
col_labels      <- c("AR=0","AR=0.5","AR=0,8","AR=0.95","AR=1")
colnames(mout1) <- col_labels
colnames(mout2) <- col_labels

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
   
    if (abs(theta - theta_old) < tol) {
      # cat("Converged in", iter, "iterations.\n")
      converged <- TRUE
      break
    }
  }
  return(list(Shared_Theta = theta, Unique_Phis = phi_list))
}

#+++++++++++++++++++++++++++++++++++
# Step 3 simulation loop (1000 reps)
#+++++++++++++++++++++++++++++++++++
phi_x <- c(0.0, 0.3, 0.5, 0.7, 0.8)
true_phi <- -0.2
for (s in 1:reps) {
  if (s %% 100 == 0) {
    cat("Now replication is", s, "out of", reps, "\n")
    flush.console()
  }
  common_y <- rnorm(T_obs)
  my_data  <- matrix(NA, nrow = T_obs, ncol = 10)
  for (j in 1:5) {
    u_x <- rnorm(T_obs)
    x <- numeric(T_obs)
    x[1] <- u_x[1]
    for (t in 2:T_obs) {
      x[t] <- phi_x[j] * x[t - 1] + u_x[t]
    }
    # Optional: demean x so it is clearly stationary around zero
    x <- x - mean(x)
    u_y <- sqrt(rho) * common_y + sqrt(1 - rho) * rnorm(T_obs)
    y <- numeric(T_obs)
    y[1] <- u_y[1]
    for (t in 2:T_obs) {
      y[t] <- y[t-1] + true_phi * (y[t-1] - true_beta * x[t]) + u_y[t]
    }
    my_data[,(j-1)*2+1]=x
    my_data[,(j-1)*2+2]=y
  }
  my_data_df <- as.data.frame(my_data)
  colnames(my_data_df) <- paste0(rep(c("x", "y"), 5), rep(1:5, each=2)) 
  my_data_df$Time_Index <- 1:T_obs
  
  panel_data <- my_data_df %>%
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
  
  results <- manual_pmg(pmg_data)
  mout1[s, ]  <- results$Unique_Phis
  mout2[s, ] <- results$Shared_Theta
  mout3[s,1]  <- round(mean(results$Unique_Phis, na.rm = TRUE), 4)
  mout4[s,1]  <- round(results$Shared_Theta, 4)
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Report the Findings
#++++++++++++++++++++++++++++++++++++++++++++++++
# Define upper-tail significance levels
alpha_levels <- c(0.01, 0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975, 0.99)

# Calculate the quantiles for all 5 equations at once
# apply(matrix, 2, ...) means "calculate this for every column"
crit_vals_alpha <- apply(mout3, 2, quantile, probs = alpha_levels, na.rm = TRUE)
crit_vals_theta <- apply(mout4, 2, quantile, probs = alpha_levels, na.rm = TRUE)

# Convert to a clean data frame
results_table_alpha <- data.frame(crit_vals_alpha)
results_table_theta <- data.frame(crit_vals_theta)


# Exactly 4 row names for the 4 significance levels
rownames(results_table_alpha) <- c("1%", "2.5%", "5%", "10%", "50%",
                             "90%", "95%", "97.5%", "99%")
rownames(results_table_theta) <- c("1%", "2.5%", "5%", "10%", "50%",
                                   "90%", "95%", "97.5%", "99%")

print("Monte Carlo Simulation Complete! GM Critical Values:")
print(round(results_table_alpha, 3))
print(round(results_table_theta, 3))
