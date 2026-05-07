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

  model <- lm(dy ~ y_lag + x_lag + dx, data = data)
  
  # Extract parameters from the model
  phi <- coef(model)["y_lag"]
  theta <- -(coef(model)["x_lag"] / phi)

  # Save the results into our storage vectors
  return(c(phi, theta))
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