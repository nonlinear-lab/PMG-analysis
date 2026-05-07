#+++++++++++++++++++++++++++++++++
# Step 0: Set working directory
#+++++++++++++++++++++++++++++++
setwd("C:/R/ardl")

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: Generate data with correlated errors
#++++++++++++++++++++++++++++++++++++++++++++++
T <- 200

# --- KEY CHANGE: shared common factor ---
# This single shock is added to ALL equations' errors.
# rho controls how strongly errors co-move across equations.
# rho = 0   --> independent errors (your original setup, mean corr ≈ 0)
# rho = 0.5 --> moderate cross-equation correlation (mean corr ≈ 0.3–0.4)
# rho = 0.8 --> high cross-equation correlation    (mean corr ≈ 0.6–0.7)

rho        <- 0.8                          # <-- change this to tune correlation
common     <- rnorm(T, mean = 0, sd = 1)   # shared shock

for (j in 1:5) {
  u_x <- rnorm(T, mean = 0, sd = 1)
  x   <- cumsum(u_x)

  # idiosyncratic (equation-specific) error component
  if (j == 1) {
    u_idio <- rnorm(T, mean = 0, sd = 1)
  } else if (j == 2) {
    u_idio <- arima.sim(model = list(ar = 0.5),  n = T, sd = 1)
  } else if (j == 3) {
    u_idio <- arima.sim(model = list(ar = 0.8),  n = T, sd = 1)
  } else if (j == 4) {
    u_idio <- arima.sim(model = list(ar = 0.95), n = T, sd = 1)
  } else {
    u_idio <- cumsum(rnorm(T, mean = 0, sd = 1))
  }

  # Combined error = rho * common shock + (1-rho) * idiosyncratic
  u_y <- rho * common + (1 - rho) * u_idio
  y   <- 10 + 2.5 * x + u_y

  assign(paste0("x", j), x)
  assign(paste0("y", j), y)
}

#+++++++++++++++++++++++++++++++++
# Step 2: Save the Data
#+++++++++++++++++++++++++++++++
master_data <- data.frame(
  Time_Index = 1:T,
  x1 = x1, y1 = y1,
  x2 = x2, y2 = y2,
  x3 = x3, y3 = y3,
  x4 = x4, y4 = y4,
  x5 = x5, y5 = y5
)
write.csv(master_data, file = "data(5).csv", row.names = FALSE)
cat("Data saved. rho =", rho, "\n")
