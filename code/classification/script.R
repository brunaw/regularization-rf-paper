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

mi_model <- function(data, gamma = 0.5, lambda_0 = 0.5){

  sets <- data %>% 
    mutate(set = ifelse(runif(n()) > 2/3, "test", "train")) %>% 
    ungroup() %>% 
    split(.$set)
  
  # Separating train and test sets in 2/3+1/3 proportion
  train <- sets$train %>% select(-set)
  test <- sets$test %>% select(-set)
  train$class <- as.factor(train$class)
  
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
             regularization = list(coef.reg = .x[[2]]))
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

# -------------
# Note: This function can still give errors and be slow. 
# -------------
# and that's it - run a 100 times each 

da_1_50 <- future_map(
  .x = 1:50, 
  ~{
    mi_model(data = final_dfs[[1]])
  })    

saveRDS(da_1_50, "rds/da_1_50.rds")


da_2_50 <- future_map(
  .x = 1:50, 
  ~{
    mi_model(data = final_dfs[[2]])
  })  

saveRDS(da_2_50, "rds/da_2_50.rds")


da_3_50 <- future_map(
  .x = 1:50, 
  ~{
    mi_model(data = final_dfs[[3]])
  })  

saveRDS(da_3_50, "rds/da_3_50.rds")

dd <- mi_model(data = final_dfs[[4]])
dd
da_4_50 <- future_map(
  .x = 1:50, 
  ~{
    mi_model(data = final_dfs[[4]])
  })  

saveRDS(da_4_50, "rds/da_4_50.rds")

da_8_50 <- future_map(
  .x = 1:50, 
  ~{
    mi_model(data = final_dfs[[8]])
  })  

saveRDS(da_8_50, "rds/da_8_50.rds")



da_3_100 <- future_map(
  .x = 1:100, 
  ~{
    mi_model(data = final_dfs[[3]])
  })  


saveRDS(da_3_100, "rds/da_3_100.rds")

plan(multiprocess)

da_4_100 <- future_map(
  .x = 1:100, 
  ~{
    mi_model(data = final_dfs[[4]])
  })

saveRDS(da_4_100, "rds/da_4_100.rds")

da_5_100 <- map(
  .x = 1:100, 
  ~{
    mi_model(data = final_dfs[[5]])
  })  

saveRDS(da_5_100, "rds/da_5_100.rds")


da_6_100 <- future_map(
  .x = 1:100, 
  ~{
    mi_model(data = final_dfs[[6]])
  })  

saveRDS(da_6_100, "rds/da_6_100.rds")




da_7_50 <- future_map(
  .x = 1:50, 
  ~{
    mi_model(data = final_dfs[[7]])
  })  

saveRDS(da_7_50, "rds/da_7_50.rds")


da_8_100 <- future_map(
  .x = 1:100, 
  ~{
    mi_model(data = final_dfs[[8]])
  })  

saveRDS(da_8_100, "rds/da_8_100.rds")

da_9_100 <- future_map(
  .x = 1:100, 
  ~{
    mi_model(data = final_dfs[[9]])
  })  

saveRDS(da_9_100, "rds/da_9_100.rds")



df <- data.frame(nvars = da_1_100 %>% 
                   map("nvars")  %>%  unlist(), 
                 error = da_1_100 %>% 
                   map("error")  %>%  unlist()) %>% 
  mutate(mtry =  rep(
    round(c(sqrt(dim(train)[2]), 
            dim(train)[2] * c(0.15, 0.40, 0.75,  0.95))), 
    5)) 

da_2_100 %>% map("mtrys")
da_2_100[[1]]
df <- data.frame(nvars = da_3_100 %>% 
                   map("nvars")  %>%  unlist(), 
                 error = da_3_100 %>% 
                   map("error")  %>%  unlist(), 
                 mtry = da_3_100 %>% 
  map("mtrys") %>% unlist())

df <- data.frame(nvars = da_1_100 %>% 
                   map("nvars")  %>%  unlist(), 
                 error = da_1_100 %>% 
                   map("error")  %>%  unlist(), 
                 mtry = da_1_100 %>% 
                   map("mtrys") %>% unlist())



df %>% 
  group_by(mtry) %>% 
  summarise(m = mean(nvars)) %>% 
  View()

df %>% 
  arrange(error)


dim(final_dfs[[4]])
df %>% 
  filter(mtry == 2290) %>% 
  group_by(nvars) %>% 
  summarise(m = median(error)) %>% 
  mutate(toff = m/nvars) %>% 
  arrange(desc(toff)) %>% 
  View()

da_1_100 %>% 
  map("error") 
  map_dbl(which.min)

  # which.max(da_1_100[[1]]$error[da_1_100[[1]]$error < 0.2]/(da_1_100[[1]]$nvars[
  #   da_1_100[[1]]$error < 0.2]))
  # 
  # model = da_1_100[[1]]
  
# group by number of parameters and find error means
  
  
find_results <- function(model){
    # ind <- which.max(model$error[model$error < 0.21]/(
    #   model$nvars[
    #     model$error < 0.21]))
    
    ords <- order(model$error)[1:2]
    ords <- c(1:5)[-(which.min(model$nvars))]
    
    ind <- which.max((model$error/model$nvars)[ords])

    error <- model$error[ords[ind]]
    nvars <- model$nvars[ords[ind]]

        # error <- model$error[model$error < 0.21][ind]
    # nvars <- model$nvars[model$error < 0.21][ind]
    # 
    # if(is_empty(ind)){
    #   ind <- which.max(model$error[model$error < 0.24]/(
    #     model$nvars[
    #       model$error < 0.24]))
    #   
    #   error <- model$error[model$error < 0.24][ind]
    #   nvars <- model$nvars[model$error < 0.24][ind]
    #   
    #   if(is_empty(ind)){
    #     ind <- which.max(model$error[model$error < 0.27]/(
    #       model$nvars[
    #         model$error < 0.27]))
    #     
    #     error <- model$error[model$error < 0.27][ind]
    #     nvars <- model$nvars[model$error < 0.27][ind]
    #   }
    #   
    # }

    return(list(error = error, nvars = nvars))  
  }  
  
da_1_100[[10]]$error
da_1_100[[10]]$nvars
res[[10]]

res %>% map("error")  %>% summary()
res <- da_1_100 %>% map(find_results)
res %>% map("error")  %>% unlist() %>% mean()
res %>% map("nvars")  %>% unlist() %>% mean()
da_1_100[[1]] %>% 
  map("nvars") %>% 
  map_dbl(min) %>% 
  mean()

data = final_dfs[[2]]


da_2_100 <- future_map(
  .x = 1:100, 
  ~{
    mi_model(data = final_dfs[[2]])
  }) 


m1 <- ranger(form, 
             data = train,
             #mtry = length(nam_sel),
             num.trees = 100)
m1

pp <- predict(m1, sets$test)
pp$predictions
1  - sum(pp$predictions == sets$test$class)/length(pp$predictions)

length(pp$predictions) * 0.08

m2 <- ranger(form, 
             data = data,
             mtry = length(nam_sel),
             num.trees = 1)
m2

models <- final_dfs %>% 
  future_map(mi_model)

length(models)

models %>% map_dbl(
  ~{.x$oob.error}
)

m1
m2

sum(m1$variable.importance > 0)
sum(m2$variable.importance > 0)

# Guided 




data = final_dfs[[1]]
m0 <- ranger(form, 
             data = final_dfs[[1]],
             num.trees = 300)


dim(final_dfs[[1]])
m1 <- ranger(form, 
             data = final_dfs[[1]],
             num.trees = 300, 
             regularization = list(coef.reg = 0.25))

sum(m0$variable.importance > 0)
sum(m1$variable.importance > 0)
m0
m1

coefs <- m0$variable.importance/max(m0$variable.importance) 

coefs_mix <- 
  cross2(gamma, lambda_0) %>% 
  map(
    ~{
      ((1 - .x[[1]]) * .x[[2]]) + .x[[1]] * coefs  
    })

map_model_guided_mix <- function(data, mtry_fc, dep_par){
  
  coefs_mix %>% 
    map(~{
      ranger(dependent.variable.name = "log_brd5", 
             data = data,
             num.trees = 200,
             regularization = list(coef_reg = .x, 
                                   use_depth  = dep_par), 
             mtry = mtry_fc, 
             num.threads = 1)
    })
}

