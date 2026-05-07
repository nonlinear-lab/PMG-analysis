#+++++++++++++++++++++++++++++++++++++++
# Step 0: Working directory and library
#+++++++++++++++++++++++++++++++++++++++
setwd("C:/R/ardl")

# Install the ardlverse package if you haven't already
# install.packages("ardlverse")

library(tidyr)      # For reshaping the data
library(dplyr)      # For data manipulation
N <- 5    # Number of groups
T <- 200  # Number of time periods

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: Load and Reshape Data
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(5).csv")

# Reshape from wide to long format
panel_data <- my_data %>%
  pivot_longer(
    cols = -Time_Index,
    names_to = c(".value", "state"),
    names_pattern = "([xy])([0-9]+)"
  )

# Convert state to numeric
panel_data$state <- as.numeric(panel_data$state)

#+++++++++++++++++++++++++++++++++++++++
# Step 1: Prepare Lags and Differences
#+++++++++++++++++++++++++++++++++++++++

# Create lagged and differenced variables required for ARDL(1,1)
# ECM format: dy_it = phi_i * (y_{i,t-1} - theta * x_{i,t}) + W_i * kappa_i + e_it
pmg_data <- panel_data %>%
  group_by(state) %>%
  mutate(
    dy    = y - lag(y),
    dx    = x - lag(x),
    y_lag = lag(y),
    x_val = x        # In the 1999 paper, the ECM is often written with current X
  ) %>%
  filter(!is.na(dy) & !is.na(dx)) %>% # Drop NA rows caused by lagging
  ungroup()

#+++++++++++++++++++++++++++++++++++++++
# Step 2: Initialize Matrices & PMG Loop
#+++++++++++++++++++++++++++++++++++++++

# Initialize the common long-run coefficient theta (can start at 0)
theta <- 0
tol <- 1e-6          # Convergence tolerance
max_iter <- 1000     # Maximum iterations
converged <- FALSE

states <- unique(pmg_data$state)

for (iter in 1:max_iter) {
  theta_old <- theta

  # Accumulators for the theta update formula
  num_sum <- 0
  den_sum <- 0

  # Lists to store group-specific estimates for the current iteration
  phi_list <- numeric(N)
  sig2_list <- numeric(N)

  # Loop over each group (state) i
  for (i in 1:length(states)) {
    group_df <- pmg_data %>% filter(state == states[i])

    # Define matrices for group i
    dy    <- as.matrix(group_df$dy)
    y_lag <- as.matrix(group_df$y_lag)
    X     <- as.matrix(group_df$x_val)

    # W contains short-run variables: differenced x and an intercept
    W     <- cbind(1, group_df$dx)

    T_i <- nrow(dy)

    # 1. Calculate the H_i Projection Matrix
    # H_i = I_T - W_i(W_i'W_i)^(-1)W_i'
    I_T <- diag(T_i)
    H_i <- I_T - W %*% solve(t(W) %*% W) %*% t(W)

    # 2. Calculate the Error Correction term (xi) given current theta
    # xi_i = y_{i,t-1} - X_i * theta
    xi <- y_lag - X * theta

    # 3. Estimate Error-Correction speed of adjustment (phi_i)
    # phi_i = (xi' H_i xi)^(-1) * (xi' H_i dy)
    phi_i <- solve(t(xi) %*% H_i %*% xi) %*% (t(xi) %*% H_i %*% dy)
    phi_i <- as.numeric(phi_i)
    phi_list[i] <- phi_i

    # 4. Estimate group-specific error variance (sigma^2_i)
    # sigma^2_i = T^(-1) * (dy - phi_i * xi)' H_i (dy - phi_i * xi)
    res_i <- dy - phi_i * xi
    sigma2_i <- as.numeric(t(res_i) %*% H_i %*% res_i) / T_i
    sig2_list[i] <- sigma2_i

    # 5. Build up the components for the new pooled theta
    # Denominator term: (phi_i^2 / sigma^2_i) * X_i' H_i X_i
    term_den <- (phi_i^2 / sigma2_i) * as.numeric(t(X) %*% H_i %*% X)
    den_sum  <- den_sum + term_den

    # Numerator term: (phi_i / sigma^2_i) * X_i' H_i (dy_i - phi_i * y_{i,t-1})
    term_num <- (phi_i / sigma2_i) * as.numeric(t(X) %*% H_i %*% (dy - phi_i * y_lag))
    num_sum  <- num_sum + term_num
  }

  # 6. Update pooled theta
  theta <- - (1 / den_sum) * num_sum

  # Check for convergence
  if (abs(theta - theta_old) < tol) {
    cat("Converged in", iter, "iterations.\n")
    converged <- TRUE
    break
  }
}

#+++++++++++++++++++++++++++++++++++++++
# Step 3: View the Results
#+++++++++++++++++++++++++++++++++++++++

if (converged) {
  cat("\n--- PMG MANUAL ESTIMATION RESULTS ---\n")
  cat("Common Long-Run Coefficient (Theta):", round(theta, 4), "\n\n")

  cat("Group-Specific Error-Correction Coefficients (Phi_i):\n")
  for (i in 1:length(states)) {
    cat("  State", states[i], ":", round(phi_list[i], 4), "\n")
  }

  # Calculate the average Error Correction speed of adjustment
  mean_phi <- mean(phi_list)
  cat("\nAverage Speed of Adjustment (Mean Phi):", round(mean_phi, 4), "\n")
} else {
  cat("Warning: Algorithm did not converge.\n")
}

