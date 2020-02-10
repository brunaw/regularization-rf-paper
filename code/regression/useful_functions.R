#-------------------------------------------------------------------
# Useful functions for the regression settings 
# 1. RMSE 
# 2. Number of parameters 
# 3. Simulation: creating data where all variables are important +
# a bunch of variables are correlated
# 4. Correlation using the Spearman method 
#-----------------------------------------------------------
library(tidyverse)
library(ranger)

# 1. Function to find the prediction power in the test set
rmse <- function(x, ntree = 100, test = test){
  
  pred <- predict(x, test, predict.all = TRUE)$predictions %>% 
    rowSums(.)/ntree
  
  sqrt(mean((pred - test$y)^2))
}

# 2. Function to calculate the number of selected variables of the model 
pars <- function(x){
  length(unique((unlist(x$forest$split.varIDs)))) 
}


# 3. Function to simulate the data as described in the paper

fried <- function(n = 1000, scale_err = 1, scale = TRUE){
  
  X = matrix(runif(n*5), nrow = n, ncol = 5)
  X_add = matrix(runif(n*200), nrow = n, ncol = 200)
  
  pars = c(0.8, 2, 1, 0.7)
  pars_add = 0.9
  
  # Creating variables with decreasing importances (0.9^(j/3))
  adds <- 1:dim(X_add)[2] %>% map(
    ~{
      (X_add[, .x ] * pars_add^(.x/3))
    }) %>% 
    bind_cols()
  
  # Creating 45 variables correlated to x_5
  adds_same <- 1:45 %>% map(
    ~{
      (X[, 5]* pars_add^(.x))
    }) %>% 
    bind_cols()
  
  # Creating the final mean
  mean = 
    (pars[1]*sin(pi* X[, 1] * X[,2]) + pars[2] * (X[,3] - 0.5)^2
     + pars[3] * X[,4] + pars[4] * X[,5]) +
    rowSums(adds)  + rowSums(adds_same)
  
  y = rnorm(n, mean, scale_err)
  
  # Saving simulated y and X and separating into train and test
  df <- data.frame(y, X) %>% 
    bind_cols(adds) %>% 
    bind_cols(adds_same) %>% 
    setNames(c("y", paste0("var_", 1:(dim(.)[2] - 1)))) %>% 
    mutate(set = ifelse(runif(nrow(.)) < 0.8, "train", "test")) %>% 
    `if`(scale == TRUE, mutate(., y = scale(y) %>% c()), .) %>% 
    split(.$set) %>% 
    map(~ .x %>% select(-set))
  
  return(list(train = df %>% pluck("train"), 
              test = df %>% pluck("test")))
  
}

# Testing it 
ff <- fried()
ff %>% map(~summary(.x$y))

# 4. Function to calculate marginal correlations in the simulated data
corr_fc <- function(data, var){
  cor(data$y, y = data %>% pull(var), method = "spearman")
}

# RRF versions  ---------------------
# 1. RMSE
rmse_rrf <- function(x, test = test){
  pred <- predict(x, test)
  sqrt(mean((pred - test$y)^2))
}
# 2. Number of parameters 
pars_rrf <- function(x){
  sum(x$importance > 0)
}

#-------------------------------------------------------------------
# End of functions
#-------------------------------------------------------------------