#-------------------------------------------------------------------
# Wrangling of the resulting models for the simulated data 
#-------------------------------------------------------------------
library(tidyverse)
library(ranger)

source("code/regression/useful_functions.R")

# Generating the parameters for all the models ---------------------
lambdas <- seq(0.05, 0.99, length.out = 15)
mtry <- c(seq(15, 250, by = 30), 250)
n <- 10

# Reading our data and separating the test set -------------------
sim_new <- readRDS("data/regression/simulated/sim_new_corr.rds")
test_list <- sim_new[2, 1:n] %>%
  rep(each = length(lambdas)*length(mtry))

# Scaled version of the data  -------------------
sim_new_sc <- sim_new[1, 1:n]  %>% 
  map(~{
    .x %>% 
      mutate(y = scale(y))
  })

test_list_sc <- sim_new[2, 1:n] %>% 
  map(~{
    .x %>% 
      mutate(y = scale(y))
  }) %>% 
  rep(each = length(lambdas)*length(mtry))

# 0. Standard Random Forests  --------------------------------
std_rf <- readRDS("rds/models/std_rf.Rds")

# 1. Constant models and scaled data --------------------------------
ct_models_scaled <- readRDS("rds/models/constant.Rds")
ct_no_depth <- ct_models_scaled$res_no_depth

#-------------------------------------------------------------------
# Building the data.frames for each case
#-------------------------------------------------------------------

# 0. Standard Random Forests --------------------

test_list_rf <- sim_new[2, 1:n] %>% 
  map(~{
    .x %>% 
      mutate(y = scale(y))
  }) %>% 
  rep(each = length(mtry))

df_rf <- data.frame(
  pars = std_rf %>% map_dbl(pars), 
  rmse =  map2_dbl(
    .x = std_rf, 
    .y = test_list_rf,
    ~{rmse(x = .x, test = .y, ntree = 100)}
  )) %>% 
  mutate(mtry_df = rep(mtry, each =  n))


# Checking importances
imps_rf <- std_rf %>% 
  map("variable.importance") %>% 
  bind_cols() %>% 
  mutate(names = names(test_list_rf[[1]][-1]))

saveRDS(list(df_rf = df_rf, 
             imps_rf = imps_rf), 
        file = "auxiliar_results/std_rf.rds")

# 1. Constant models --------------------
df_ct <- 
  # Without consider the depth 
  data.frame(
    pars = c(map_dbl(ct_no_depth, pars)),
    rmse = 
      map2_dbl(
        .x = ct_no_depth, 
        .y = test_list_sc,
        ~{rmse(x = .x, test = .y, ntree = 100)}
      )) %>% 
  mutate(type_depth = rep(c("No Depth"),
                          each = length(lambdas)*length(mtry)*n), 
         lambdas_df = rep(rep(lambdas, each = n), times = length(mtry)), 
         mtry_df  = rep(mtry, each = length(lambdas) * n))

write.table(df_ct, file = "auxiliar_results/df_ct.txt")

rm(ct_no_depth)
rm(ct_models_scaled)
rm(std_rf)

#-------------------------------------------------------------------
# Mixture models
#-------------------------------------------------------------------
gamma <- seq(0.001, 0.999, length.out = 5)
lambda_0 <- seq(0.1, 0.9, length.out = 5)

test_list_sc <- sim_new[2, 1:n] %>% 
  map(~{
    .x %>% 
      mutate(y = scale(y))
  }) %>% 
  rep(each = length(mtry))


# Reading models back: 
boosted_rf <- readRDS("rds/models/boosted_rf.Rds")
corr_models <- readRDS("rds/models/corr_models.Rds")
boosted_svm <- readRDS("rds/models/boosted_svm.Rds")

# 3. Boosted RF
df_boosted_rf <- 
  # Without consider the depth 
  data.frame(
    pars = boosted_rf$boosted_rf %>% 
      map(
        ~{ .x %>% map_dbl(pars)}) %>% unlist(),
    rmse = 
      boosted_rf$boosted_rf %>% 
      map2(
        .y = test_list_sc, 
        ~{.x %>% map_dbl(rmse, test = .y)})  %>% 
      unlist() 
  ) %>% 
  mutate(
    
    lambdas_df = rep(rep(lambda_0, each = length(gamma)),
                     times = length(mtry)*10),
    gamma_df = rep(gamma, length(mtry) * 10 * length(lambda_0)), 
    mtry_df  = rep(
      rep(mtry, each = length(lambda_0) *  length(gamma)), 10) 
  ) %>% 
  mutate(dataset = rep(1:10, each = length(lambda_0) * length(gamma)* length(mtry))) %>% 
  mutate(model = "boosted_rf")


# 2. Correlation
df_corr <- 
  # Without consider the depth 
  data.frame(
    pars = corr_models$corr_models %>% 
      map(
        ~{ .x %>% map_dbl(pars)}) %>% unlist(),
    rmse = 
      corr_models$corr_models %>% 
      map2(
        .y = test_list_sc, 
        ~{.x %>% map_dbl(rmse, test = .y)})  %>% 
      unlist() 
  ) %>% 
  mutate(
    
    lambdas_df = rep(rep(lambda_0, each = length(gamma)),
                     times = length(mtry)*10),
    gamma_df = rep(gamma, length(mtry) * 10 * length(lambda_0)), 
    mtry_df  = rep(
      rep(mtry, each = length(lambda_0) *  length(gamma)), 10) 
  ) %>% 
  mutate(dataset = rep(1:10, each = length(lambda_0) * length(gamma)* length(mtry))) %>% 
  mutate(model = "correlation")


# 3. Boosted SVM

df_boosted_svm <- 
  # Without consider the depth 
  data.frame(
    pars = boosted_svm$boosted_svm %>% 
      map(
        ~{ .x %>% map_dbl(pars)}) %>% unlist(),
    rmse = 
      boosted_svm$boosted_svm %>% 
      map2(
        .y = test_list_sc, 
        ~{.x %>% map_dbl(rmse, test = .y)})  %>% 
      unlist() 
  ) %>% 
  mutate(
    
    lambdas_df = rep(rep(lambda_0, each = length(gamma)),
                     times = length(mtry)*10),
    gamma_df = rep(gamma, length(mtry) * 10 * length(lambda_0)), 
    mtry_df  = rep(
      rep(mtry, each = length(lambda_0) *  length(gamma)), 10) 
  ) %>% 
  mutate(dataset = rep(1:10, 
                       each = length(lambda_0) * length(gamma)* length(mtry))) %>% 
  mutate(model = "boosted_svm")


df_mix <- df_boosted_rf  %>% 
  bind_rows(df_corr) %>% 
  bind_rows(df_boosted_svm)


# Extracting variable names used ---
nam <- names(sim_new_sc[[1]])

get_names <- function(model){
  index <- unique(unlist(model$forest$split.varIDs))
  nam[index]
}

nam_rf <- boosted_rf$boosted_rf %>% 
  map( ~{.x %>% map(get_names)} ) 

nam_corr <- corr_models$corr_models %>% 
  map(~{ .x %>% map(get_names)}) 

nam_svm <- boosted_svm$boosted_svm %>% 
  map(~{ .x %>% map(get_names)}) 

# Finding proportions of important and of correlated variables ---
vec <- 0.9^(1:200/3)

limi <- sum(vec > 0.01)
imp_nam <- nam[-1][c(1:5, 6:limi)]
corr_vars <- nam[-1][c(206:250)]
vec <- 0.9^(1:45)

# Identifying proportions
identify_imp <- function(names_in_model){
  
  using_corr <- sum(corr_vars %in% names_in_model)
  using_imp <- sum(imp_nam %in% names_in_model[
    which(!names_in_model %in% corr_vars)])
  
  if(using_imp > 0){
    good <- using_imp/(length(names_in_model) - using_corr) }
  else{ good <- 0 }
  
  bad <- using_corr/length(corr_vars)
  
  return(list(good = good, bad = bad))
}

imps_rf <- nam_rf %>% 
  map(~{.x %>% map(identify_imp)}) 

imps_corr <- nam_corr %>% 
  map(~{ .x %>% map(identify_imp)}) 

imps_svm <- nam_svm %>% 
  map(~{ .x %>% map(identify_imp)}) 

# Saving results --- 
df_mix$good <- c(imps_rf, imps_corr, imps_svm) %>%
  map(~map_dbl(.x, "good")) %>% unlist()

df_mix$bad <- c(imps_rf, imps_corr, imps_svm) %>%
  map(~map_dbl(.x, "bad")) %>% unlist()

# Saving it --- 
write.table(df_mix, file = "auxiliar_results/df_mix.txt")

rm(boosted_rf)
rm(corr_models)
rm(boosted_svm)

#-------------------------------------------------------------------
# End of wrangling for results with the ranger package
#-------------------------------------------------------------------

