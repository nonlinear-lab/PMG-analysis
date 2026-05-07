#+++++++++++++++++++++++++++++++++++++++
# Step 0: dorking directory and library
#+++++++++++++++++++++++++++++++++++++++
setwd("C:/R/ardl")
# install.packages("plm")
# Load the library
library(plm)
library(tidyr) # For reshaping the data
library(dplyr) # For data manipulation

mout    <- matrix(NA, nrow=2, ncol=1)

#++++++++++++++++++++++++++++++++++++++++++++++++
# step1; estimation static model
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(5).csv")

panel_data <- my_data %>%
  pivot_longer(
    cols = -Time_Index,                  # Pivot all columns except Time_Index
    names_to = c(".value", "state"),     # Split the column names into variable (x/y) and state ID (1-5)
    names_pattern = "([xy])([0-9]+)"     # Regex: letters 'x' or 'y' go to .value, numbers go to 'state'
  )
# Convert state to a factor/numeric if needed
panel_data$state <- as.numeric(panel_data$state)

# 1. Mean Group (MG) Estimation
# We specify model = "mg" to calculate the Mean Group estimator.
# The index argument defines the group and time dimensions.
mg_model <- pmg(y ~ x, data = panel_data,
                index = c("state", "Time_Index"),
                model = "mg")
cat("\n======================================================\n")
cat("   RESULTS 1: STATIC MEAN GROUP MODEL (No Dynamics)   \n")
cat("======================================================\n")
print(summary(mg_model))

#++++++++++++++++++++++++++++++++++++++++++++++++
# step1; estimation static model
#++++++++++++++++++++++++++++++++++++++++++++++++

# First, create a pdata.frame so R understands the panel lags/differences
p_data <- pdata.frame(panel_data, index = c("state", "Time_Index"))

# Run the pmg using the Error Correction formula!
mg_ecm <- pmg(diff(y) ~ lag(y) + lag(x) + diff(x),
              data = p_data,
              model = "mg")

cat("\n======================================================\n")
cat("  RESULTS 2: DYNAMIC ECM MEAN GROUP MODEL (With Phi)  \n")
cat("======================================================\n")
print(summary(mg_ecm))

mout[1,1] <- coef(mg_model)["x"]
mout[2,1] <- coef(mg_ecm)["lag(y)"]

print(mout)
