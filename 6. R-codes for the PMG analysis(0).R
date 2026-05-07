#+++++++++++++++++++++++++++++++++++++++
# Step 0: Working directory and library
#+++++++++++++++++++++++++++++++++++++++
setwd("C:/R/ardl")

# Install the ardlverse package if you haven't already
# install.packages("ardlverse")

library(tidyr)      # For reshaping the data
library(dplyr)      # For data manipulation
library(ardlverse)  # For Panel ARDL / PMG estimation

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

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Pooled Mean Group (PMG) Estimation
#++++++++++++++++++++++++++++++++++++++++++++++++
# Using the panel_ardl() function, we set estimator = "pmg".
# We must also specify the number of lags for y (p) and x (q). 
# We'll use p=1, q=1 for a standard ARDL(1,1) model.

pmg_model <- panel_ardl(
  formula = y ~ x, 
  data = panel_data,
  id = "state",           # The column identifying the groups
  time = "Time_Index",    # The column identifying the time periods
  p = 1,                  # Lag order for the dependent variable (y)
  q = 1,                  # Lag order for the independent variable(s) (x)
  estimator = "pmg"       # Specify Pooled Mean Group
)

# View the PMG coefficient estimates (Long-run, Short-run, and Speed of Adjustment)
summary(pmg_model)