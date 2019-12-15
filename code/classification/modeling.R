library(ranger)
library(stringr)
library(tidyverse)
library(furrr)
library(glmnet)
library(varSelRF)
library(infotheo) # for MI package
plan(multiprocess)

# Data reading ---------------------------------
# Gene datasets 

ds <- list.files("data/classification/varselif", 
                 pattern = "data\\.txt$") %>% 
  .[-7] %>% 
  map(~{
    file <- paste0("data/classification/varselif/", .x, collapse = "")
    d <- read_table2(file) 
    if("#" %in% names(d)){
      d <- d %>% select(-1)
    }
    d <- d %>% 
      data.table::transpose() %>% 
      as.data.frame()
    })

# Putting back the classes
classes <- list.files("data/classification/varselif", 
                 pattern = "class\\.txt$") %>% 
  map(~{
    file <- paste0("data/classification/varselif/", .x, collapse = "")
    scan(file, character(), quote = "") 
  })

# Saving final data.frames in an object
final_dfs <- map2(
  .x = ds, .y = classes, 
  ~{
    .x %>% 
      mutate(class = .y)
  })

data = final_dfs[[2]]
# Models ---------------------------------
# 1. Boosted MI
# 2. Boosted RF
# 3. Standard RF


# Settings ------
# MI function
mutif_fc <- function(data, var){
  mutinformation(data$class %>% c(), data %>% pull(var))
}

# Parameters: gamma = 0.5; lambda_0 = 0.5 
# (half-half mixture)

mi_model <- function(train, test, gamma = 0.5, lambda_0 = 0.5){
  
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
  # 1. Boosted MI
  coefs_mix[[1]] <- ((1 - gamma) * lambda_0 + 
                       gamma * (mi_vars/max(mi_vars)))
  
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


  # Defining mtry values: from sqrt(p) to 95% of the variables 
  props <- c(round(sqrt(dim(train)[2])/dim(train)[2], 2), 0.15, 0.40, 0.75,  0.95)
  mtrys <- round(props * dim(train)[2])
  
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
  
 # varSelRF 
  xdata = train %>% select(-class) %>% as.matrix() 
  
  varselrf <- props %>% 
    map(~{
      varSelRF(xdata = xdata, 
               Class = factor(train$class), 
               ntree = 500, 
               mtryFactor = .x)
    })
  
  # LASSO
  glmmod <- cv.glmnet(xdata, y = as.factor(train$class),
                      alpha = 1, family = "multinomial", 
                      type.measure="auc")
  
  # Saving all models --------------
  mod <- models %>% prepend(varselrf)
  mod2 <- mod %>% prepend(list(glmmod))

  # Formatting results 
  form_models <- mod2 %>% 
    map(~{
      if(class(.x) == "varSelRF"){
        nam_sel <- .x$selected.vars
      } 
      else if(class(.x) == "cv.glmnet"){
        cc <- coef(.x, s = "lambda.min")
        ind <- cc %>% map( ~{which(.x != 0)}) %>% unlist() %>% unique() 
        
        nam_sel <- names(train)[-dim(train)[2]][ind]
        
      } else {
        nam_sel <- names(train)[-dim(train)[2]][
          .x$variable.importance > 0
          ]}
      
      
      form <- paste0('class ~ ', 
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
              model = c(
                "LASSO", rep(c("varSelRF", "MI", "GRRF", "RF"), 
                             each = length(mtrys)))
  ))
}


smi_model <- safely(mi_model)

# -----------------------
# Note: This function can still give errors and be slow. 

# Had to do separations by hand because of weird errors (still)

# Separating train and test sets in 2/3+1/3 proportion

# -------------
# 1 = too many problems; 
# data <- final_dfs[[1]]
# 
# sets <- data %>% 
#   mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
#   ungroup() %>% 
#   split(.$set)
# 
# train <- sets$train %>% select(-set)
# test <- sets$test %>% select(-set)
# train$class <- as.factor(train$class)
# 
# da_1_10 <- future_map(
#   .x = 1:10, 
#   ~{
#     mi_model(train = train, test = test)
#   })    
# 
# saveRDS(da_1_10, "rds/da_1_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_2_10, "rds/da_2_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_3_10, "rds/da_3_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_4_10, "rds/da_4_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_6_10, "rds/da_6_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_7_10, "rds/da_7_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_8_10, "rds/da_8_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_9_10, "rds/da_9_10.rds")
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
    mi_model(train = train, test = test)
  })    

saveRDS(da_10_10, "rds/da_10_10.rds")
# --------------------------------------------------------
