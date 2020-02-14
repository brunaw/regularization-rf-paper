library(ranger)
library(tidyverse)
library(furrr)
library(infotheo) # for MI package
# Parallelization (!!!)
plan(multiprocess) 

# Data reading ---------------------------------
# Gene datasets 

ds <- list.files("data/classification/files", 
                 pattern = "data\\.txt$") %>% 
  .[-7] %>% 
  map(~{
    file <- paste0("data/classification/files/", .x, collapse = "")
    d <- read_table2(file) 
    if("#" %in% names(d)){
      d <- d %>% select(-1)
    }
    d <- d %>% 
      data.table::transpose() %>% 
      as.data.frame()
  })

# Putting back the classes
classes <- list.files("data/classification/files", 
                      pattern = "class\\.txt$") %>% 
  map(~{
    file <- paste0("data/classification/files/", .x, collapse = "")
    scan(file, character(), quote = "") 
  })

# Saving all final data.frames in an object
final_dfs <- map2(
  .x = ds, .y = classes, 
  ~{
    .x %>% 
      mutate(class = .y)
  })


final_dfs[[1]] <-  final_dfs[[1]][-1,]  # This row is all characters 

# Saving and reading it back
final_dfs %>% saveRDS("data/classification/final_dfs.Rds")

final_dfs <-  readRDS("data/classification/final_dfs.Rds")


# Models ---------------------------------
# 0. GRRF with the RRF package
# 1. Boosted MI
# 2. Boosted RF
# 3. Standard RF

# Settings ------
# Mutual Information calculation function
mutif_fc <- function(data, var){
  mutinformation(data$class %>% c(), data %>% pull(var))
}

# Parameters: gamma = 0.5; lambda_0 = 0.5 
# (half-half mixture) for all the different g(x_i)

# Function without internal sampling
# Not used in the paper!

modelling_classification <- function(
  train, test, gamma = 0.5, lambda_0 = 0.5){
  
  # Defining mtry values: from sqrt(p) to 95% of the variables 
  props <- c(round(sqrt(dim(train)[2])/dim(train)[2], 2), 0.15, 0.40, 0.75,  0.95)
  mtrys <- round(props * dim(train)[2])
  
  # Calculating the mutual information values
  disc_data  <- discretize(train) 
  disc_data$class <- as.factor(train$class)
  nam <- names(train)[-dim(train)[2]]
  
  mi_vars <- nam  %>% map_dbl(~{ mutif_fc(data = disc_data, var = .x) })
  
  # Setting the formula
  form <- class ~ .
  
  # Calculating the equivalent penalization coefficients
  coefs_mix <- list()
  
  # 0. GRRF
  rf_model <- RRF::RRF(form, data = train, 
                       ntree = 500, 
                       flagReg = 0)
  imp_grrf <- (
    (1 - gamma) + gamma * (rf_model$importance/max(rf_model$importance)))
  
  # Needs to be mapped by mtrys
  grrf_models <-  mtrys %>% 
    map(~{
      RRF::RRF(form, data = train, 
               mtry = .x,
               coefReg = imp_grrf, 
               flagReg = 1)
    })
  
  # 1. Boosted MI
  coefs_mix[[1]] <- (
    (1 - gamma) * lambda_0 + gamma * (mi_vars/max(mi_vars)))
  
  # 2. Boosted RF
  mm <-  ranger(form, 
                data = train,
                num.trees = 500, 
                importance = "impurity")
  
  coefs_mix[[2]] <- 
    ((1 - gamma) * lambda_0 + gamma * mm$variable.importance/max(mm$variable.importance))
  
  # Standard RF
  coefs_mix[[3]] <- 1
  
  # Constraints
  coefs_mix <- coefs_mix %>% 
    map(~{
      .x %>% 
        ifelse(. > 1, 1, .) %>% 
        ifelse(. < 0, 0.00001, .) 
    })
  
  # Combining all parameters and running the models
  models <- mtrys %>% 
    cross2(coefs_mix) %>% 
    map(~{
      ranger(form, 
             data = train,
             mtry = .x[[1]], 
             num.trees = 500, 
             regularization.factor = .x[[2]])
    })
  
  # Saving all models --------------
  models <- models %>% prepend(grrf_models)
  
  # Formatting results 
  form_models <- models %>% 
    map(~{
      
      # For GRRF
      if("RRF" %in% class(.x)){
        nam_sel <- names(train)[-dim(train)[2]][.x$importance > 0]
        
      } else{
        
        # For the models using the ranger function
        nam_sel <- names(train)[-dim(train)[2]][
          .x$forest$split.varIDs %>% unlist() %>% unique() + 1
          ]
        }
      
      form <- paste0("class ~ ", paste0(nam_sel, collapse = '+')) %>% 
        as.formula()
      
      m1 <- ranger(form, 
                   data = train,
                   num.trees = 100)
      
      pp <- predict(m1, sets$test)
      error <- 1  - sum(pp$predictions == test$class)/length(pp$predictions)
      list(error = error, nvars = length(nam_sel))
    })
  
  # Returning everything -----
  return(list(error = form_models %>% map_dbl("error"), 
              nvars = form_models %>% map_dbl("nvars"), 
              mtrys = mtrys, 
              model = rep(c("Boosted_MI", "Boosted_RF", "Std_RF", "GRRF"), 
                          each = length(mtrys)))
  )
}

# This should prevent the function to stop running if an error 
# is found
modelling_classification <- safely(modelling_classification)

# Separating train and test sets in 2/3+1/3 proportion -------
data <- final_dfs[[1]]

sets <- data %>%
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>%
  ungroup() %>%
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_1_10 <- future_map(
  .x = 1:10,
  ~{ modelling_classification(train = train, test = test) })

saveRDS(da_1_10, "auxiliar_results/classification/da_1_10.rds")
# -------------
# 2 dataset
data <- final_dfs[[2]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_2_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_2_10, "auxiliar_results/classification/da_2_10.rds")
# -------------
# 3 dataset
data <- final_dfs[[3]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)


da_3_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_3_10, "auxiliar_results/classification/da_3_10.rds")
# -----------------------
# 4 dataset
data <- final_dfs[[4]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_4_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_4_10, "auxiliar_results/classification/da_4_10.rds")
# -----------------------
# 4 dataset
data <- final_dfs[[5]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_5_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_5_10, "auxiliar_results/classification/da_5_10.rds")
# -----------------------
# 6 dataset
data <- final_dfs[[6]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_6_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_6_10, "auxiliar_results/classification/da_6_10.rds")
# -----------------------
# 7 dataset
data <- final_dfs[[7]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_7_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_7_10, "auxiliar_results/classification/da_7_10.rds")
# -----------------------
# 8 dataset
data <- final_dfs[[8]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_8_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_8_10, "auxiliar_results/classification/da_8_10.rds")
# -----------------------
# 9 dataset
data <- final_dfs[[9]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_9_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_9_10, "auxiliar_results/classification/da_9_10.rds")
# -----------------------
# 10 dataset
data <- final_dfs[[10]]

sets <- data %>% 
  mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
  ungroup() %>% 
  split(.$set)

train <- sets$train %>% select(-set)
test <- sets$test %>% select(-set)
train$class <- as.factor(train$class)

da_10_10 <- future_map(
  .x = 1:10, 
  ~{
    modelling_classification(train = train, test = test)
  })    

saveRDS(da_10_10, "auxiliar_results/classification/da_10_10.rds")

# --------------------------------------------------------
# With internal sampling
# This is what is used in the paper!
# --------------------------------------------------------

modelling_classification_with_sampling <- function(
  data, gamma = 0.5, lambda_0 = 0.5){
  
  # Now the sampling is being done inside the function, for each 
  # time it is run
  sets <- data %>% 
    group_by(class) %>% 
    mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
    ungroup() %>% 
    split(.$set)
  
  train <- sets$train %>% select(-set)
  test <- sets$test %>% select(-set)
  train$class <- as.factor(train$class)
  
  # Defining mtry values: from sqrt(p) to 95% of the variables 
  props <- c(round(sqrt(dim(train)[2])/dim(train)[2], 2), 0.15, 0.40, 0.75,  0.95)
  mtrys <- round(props * dim(train)[2])
  
  # Calculating the mutual information values
  disc_data  <- discretize(train) 
  disc_data$class <- as.factor(train$class)
  nam <- names(train)[-dim(train)[2]]
  
  mi_vars <- nam  %>% 
    map_dbl(~{ mutif_fc(data = disc_data, var = .x) })
  
  # Setting the formula
  form <- class ~ .
  
  # Calculating the equivalent penalization coefficients
  coefs_mix <- list()
  
  
  # 0. GRRF
  rf_model <- RRF::RRF(form, data = train, 
                       ntree = 500, 
                       flagReg = 0)
  imp_grrf <- (
    (1 - gamma) +  gamma * (rf_model$importance/max(rf_model$importance)))
  
  # Needs to be mapped by mtry
  grrf_models <-  mtrys %>% 
    map(~{
      RRF::RRF(form, data = train, 
               mtry = .x,
               coefReg = imp_grrf, 
               flagReg = 1)
    })
  
  # 1. Boosted MI
  coefs_mix[[1]] <- (
    (1 - gamma) * lambda_0 + gamma * (mi_vars/max(mi_vars)))
  
  # 2. Boosted RF
  mm <-  ranger(form, 
                data = train,
                num.trees = 500, 
                importance = "impurity")
  
  coefs_mix[[2]] <- 
    ((1 - gamma) * lambda_0 + gamma * mm$variable.importance/max(mm$variable.importance))
  
  # Standard RF
  coefs_mix[[3]] <- 1
  
  # Constraints
  coefs_mix <- coefs_mix %>% 
    map(~{
      .x %>% 
        ifelse(. > 1, 1, .) %>% 
        ifelse(. < 0, 0.00001, .) 
    })
  
  
  # Crossing all parameters and modelling
  models <- mtrys %>% 
    cross2(coefs_mix) %>% 
    map(~{
      ranger(form, 
             data = train,
             mtry = .x[[1]], 
             num.trees = 500, 
             regularization.factor = .x[[2]])
    })
  
  # Saving all models --------------
  models <- models %>% prepend(grrf_models)
  
  # Formatting results 
  form_models <- models %>% 
    map(~{
      # For GRRF 
      if("RRF" %in% class(.x)){
        nam_sel <- names(train)[-dim(train)[2]][.x$importance > 0]
        
      } else{
        # For ranger
        nam_sel <- names(train)[-dim(train)[2]][
          .x$forest$split.varIDs %>% unlist() %>% unique() + 1
          ]
      }
      
      form <- paste0("class ~ ", 
                     paste0(nam_sel, collapse = '+')) %>% 
        as.formula()
      
      m1 <- ranger(form, 
                   data = train,
                   num.trees = 100)
      
      pp <- predict(m1, sets$test)
      error <- 1  - sum(pp$predictions == test$class)/length(pp$predictions)
      list(error = error, nvars = length(nam_sel))
    })
  
  # Returning everything -----
  return(list(error = form_models %>% map_dbl("error"), 
              nvars = form_models %>% map_dbl("nvars"), 
              mtrys = mtrys, 
              model = rep(c("Boosted_MI", "Boosted_RF", "Std_RF", "GRRF"), 
                          each = length(mtrys)))
  )
}

modelling_classification <- safely(modelling_classification)
# --------------------------------------------------------
# Running all the models again, 50 times
da_1_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[1]]) })

saveRDS(da_1_ws, "auxiliar_results/classification/da_1_ws.rds")
# --------------------------------------------------------
da_2_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[2]]) })

saveRDS(da_2_ws, "auxiliar_results/classification/da_2_ws.rds")
# --------------------------------------------------------
da_3_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[3]]) })

saveRDS(da_3_ws, "auxiliar_results/classification/da_3_ws.rds")
# --------------------------------------------------------
da_4_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[4]]) })

saveRDS(da_4_ws, "auxiliar_results/classification/da_4_ws.rds")
# --------------------------------------------------------
da_5_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[5]]) })

saveRDS(da_5_ws, "auxiliar_results/classification/da_5_ws.rds")
# --------------------------------------------------------
da_6_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[6]]) })

saveRDS(da_6_ws, "auxiliar_results/classification/da_6_ws.rds")
# --------------------------------------------------------
da_7_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[7]]) })

saveRDS(da_7_ws, "auxiliar_results/classification/da_7_ws.rds")
# --------------------------------------------------------
da_8_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[8]]) })

saveRDS(da_8_ws, "auxiliar_results/classification/da_8_ws.rds")
# --------------------------------------------------------
da_9_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[9]]) })

saveRDS(da_9_ws, "auxiliar_results/classification/da_9_ws.rds")
# --------------------------------------------------------
da_10_ws <- future_map(
  .x = 1:50,
  ~{ modelling_classification_with_sampling(data = final_dfs[[10]]) })

saveRDS(da_10_ws, "auxiliar_results/classification/da_10_ws.rds")
# -------------------------------------------------------------