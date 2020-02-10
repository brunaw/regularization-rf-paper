#----------------------------------------------------------
# Mixture models for the simulated data 
#----------------------------------------------------------
library(tidyverse)
library(ranger)
library(furrr)
library(rminer)
plan(multiprocess)

source("code/regression/useful_functions.R")

# Generating the parameters for all the models ---------------------
lambdas <- seq(0.05, 0.99, length.out = 15)
gammas_grrf <- lambdas
mtry <- c(seq(15, 250, by = 30), 250)
n <- 10
gamma <- seq(0.001, 0.999, length.out = 5)
lambda_0 <- seq(0.1, 0.9, length.out = 5)

# Generating our data and separating the test set -------------------
# Run only once
# sim_new <- replicate(n = n, fried(scale = FALSE)) 
# write_rds(sim_new, "data/regression/simulated/sim_new_corr.rds")
sim_new <- readRDS("data/regression/simulated/sim_new_corr.rds")
test_list <- sim_new[2, 1:n] %>%
  rep(each = length(lambdas)*length(mtry))

# Scaled version of the data  -------------------
sim_new_sc <- sim_new[1, 1:n]  %>% 
  map(~{
    .x %>% mutate(y = scale(y))
  })

test_list_sc <- sim_new[2, 1:n] %>% 
  map(~{
    .x %>% mutate(y = scale(y))
  }) %>% 
  rep(each = length(lambdas)*length(mtry))

#-------------------------------------------------------------------
# 0. Standard random forest
#-------------------------------------------------------------------

std_rf <- sim_new_sc[1:n] %>% 
  purrr::cross2(mtry) %>% 
  furrr::future_map(
    ~{
      ranger(y ~ ., 
             data = .x[[1]], 
             mtry = .x[[2]], 
             num.trees = 100,
             importance = "impurity")
    })

saveRDS(std_rf, file = "rds/models/std_rf.Rds")
#-------------------------------------------------------------------
# 1. Constant \lambda_i models
#-------------------------------------------------------------------

# Building the modelling function 
fc_mod <- function(data, coef_reg, dep_par,
                   mtry, num.trees = 100){
  ranger(y ~ ., 
         data = data,
         num.trees = num.trees,
         mtry = mtry,
         importance = "impurity",
         regularization.factor  = coef_reg,
         regularization.usedepth = dep_par)
}

# Without the depth 
res_no_depth <- sim_new_sc[1:n] %>% 
  purrr::cross3(lambdas, mtry) %>%
  furrr::future_map(~fc_mod(data = .x[[1]], 
                            coef_reg = .x[[2]],
                            mtry = .x[[3]],
                            dep_par = FALSE))

saveRDS(list(sim_sc = sim_new_sc, 
             res_no_depth = res_no_depth), 
        file = "rds/models/constant.Rds")
#-------------------------------------------------------------------
# 2.1 Constant \lambda_i models with the RRF package
#-------------------------------------------------------------------

# Building the modelling function 
fc_mod_rrf <- function(data, coefReg, 
                       mtry, num.trees = 100){
  RRF::RRF(y ~ ., 
           data = data,
           coefReg  = coefReg,
           flagReg = 1, 
           ntree = num.trees,
           mtry = mtry)
}

res_rrf <- sim_new_sc[1:n] %>% 
  purrr::cross3(lambdas, mtry) %>%
  purrr::map(~fc_mod_rrf(data = .x[[1]], 
                         coefReg = .x[[2]],
                         mtry = .x[[3]]))

df_rrf_scaled <-
  data.frame(
    pars = c(map_dbl(res_rrf, pars_rrf)),
    rmse = 
      map2_dbl(
        .x = res_rrf, 
        .y = test_list_sc,
        ~{rmse_rrf(x = .x, test = .y)}
      )) %>% 
  mutate(type_depth = rep(c("No Depth"),each = length(lambdas)*length(mtry)*n), 
         lambdas_df = rep(rep(lambdas, each = n), times = length(mtry)), 
         mtry_df  = rep(mtry, each = length(lambdas) * n)
  )

saveRDS(list(df_rrf_scaled = df_rrf_scaled), 
        file = "auxiliar_results/df_rrf_scaled_new.Rds")
#-------------------------------------------------------------------
# 2.2 GRRF models with the RRF package
#-------------------------------------------------------------------
test_list_grrf <- sim_new[2, 1:n] %>% 
  map(~{
    .x %>% mutate(y = scale(y))
  }) %>% 
  rep(each = length(mtry))

fc_mod_grrf <- function(data, mtry, num.trees = 100){
  m0 <- RRF::RRF(y ~ ., data = data, flagReg = 0)
  imps <- m0$importance/max(m0$importance)
  
  coefs_mix <- 
    gammas_grrf %>% 
    map(~{ (1 - .x) + .x * imps })
  
  coefs_mix %>% 
    map(~{
      RRF::RRF(y ~ ., 
               data = data, 
               ntrees = 100, 
               mtry = mtry, 
               coefReg = .x, 
               flagReg = 1)
    })
}
  
  
res_grrf <- sim_new_sc[1:n] %>% 
  purrr::cross2(mtry) %>%
  furrr::future_map(~fc_mod_grrf(data = .x[[1]], mtry = .x[[2]]))

df_grrf_scaled <-
  data.frame(
    pars = map(res_grrf, ~{.x %>% map_dbl(pars_rrf)}) %>% unlist(),
    rmse = 
        res_grrf  %>% 
          map2(
            .y = test_list_grrf, 
            ~{.x %>% map_dbl(rmse_rrf, test = .y)})  %>% 
          unlist() 
      ) %>% 
  mutate( 
         lambdas_df = rep(gammas_grrf,  times = length(mtry)*n), 
         mtry_df  = rep(mtry, each = length(gammas_grrf) * n))


find_vars_rrf <- function(model){
  names(test_list_grrf[[1]])[-1][model$importance > 0]  
}

vars_grrf <-  map(res_grrf$res_grrf, ~{.x %>% map(find_vars_rrf)}) 

saveRDS(list(df_grrf_scaled = df_grrf_scaled), 
        file = "auxiliar_results/df_grrf_scaled_new.Rds")

saveRDS(list(vars_grrf = vars_grrf), 
        file = "auxiliar_results/vars_grrf.Rds")

saveRDS(list(res_grrf = res_grrf), file = "rds/res_grrf.Rds")

#-------------------------------------------------------------------
# 3. Boosted RF
#-------------------------------------------------------------------

# Function to map data over the models ----------------------------
map_boosted_rf <- function(data, mtry_fc, dep_par){
  
  m0 <- ranger(y ~ ., 
               data = data, importance = "impurity")
  coefs <- m0$variable.importance/max(m0$variable.importance) 
  
  coefs_mix <- 
    cross2(gamma, lambda_0) %>% 
    map(~{ ((1 - .x[[1]]) * .x[[2]]) + .x[[1]] * coefs  })
  
  # Mapping for all the combinations of parameters
  coefs_mix %>% 
    map(~{
      ranger(y ~ .,
             data = data, 
             mtry = mtry_fc, 
             importance = "impurity",
             num.trees = 100,
             num.threads = 1,  # no parallelization 
             regularization.factor = .x, 
             regularization.usedepth = dep_par)
    })
}

boosted_rf <- mtry %>% 
  cross2(sim_new_sc) %>% 
  future_map(
    ~{
      map_boosted_rf(
        data = .x[[2]],
        mtry_fc = .x[[1]], 
        dep_par = FALSE
      )
    })

list(boosted_rf = boosted_rf) %>% saveRDS("rds/models/boosted_rf.Rds")

#-------------------------------------------------------------------
# 4. Correlation
#-------------------------------------------------------------------

# Calculating all the marginal correlations in all 10 datasets
corr_fc <- Vectorize(corr_fc, "var")

map_corr <- function(data, mtry_fc, dep_par){
  
  nam <- names(data)[-1]
  corr_vars <- nam %>% 
    map_dbl(~{corr_fc(data = data, var = .x) %>% abs()})
  
  coefs_mix <- 
    cross2(gamma, lambda_0) %>% 
    map(~{ ((1 - .x[[1]]) * .x[[2]]) + .x[[1]] * corr_vars[[1]]  })
  
  # Mapping for all the combinations of parameters
  coefs_mix %>% 
    map(~{
      ranger(y ~ .,
             data = data, 
             mtry = mtry_fc, 
             importance = "impurity",
             num.trees = 100,
             num.threads = 1,  # no parallelization 
             regularization.factor = .x, 
             regularization.usedepth = dep_par)
    })
}

corr_models <- mtry %>% 
  cross2(sim_new_sc) %>% 
  future_map(
    ~{
      map_corr(
        data = .x[[2]],
        mtry_fc = .x[[1]], 
        dep_par = FALSE
      )})

list(corr_models = corr_models) %>% saveRDS("rds/models/corr_models.Rds")

#-------------------------------------------------------------------
# 5. Boosted SVM 
#-------------------------------------------------------------------

map_boosted_svm <- function(data, mtry_fc, dep_par){
  
  M <- fit(y ~ ., data=data, 
           model="svm", kpar=list(sigma=0.10),  C = 2)
  svm.imp <- Importance(M, data=data)
  imp_norm <- (svm.imp$imp[-1]/max(svm.imp$imp[-1]))
  
  check_corr <- imp_norm == lag(imp_norm)
  check_corr[is.na(check_corr)] <- TRUE
  imp_norm[check_corr] <- 0
  
  coefs_mix <- 
    cross2(gamma, lambda_0) %>% 
    map(~{((1 - .x[[1]]) * .x[[2]]) + .x[[1]] * imp_norm})
  
  # Mapping for all the combinations of parameters
  coefs_mix %>% 
    map(~{
      ranger(y ~ .,
             data = data, 
             mtry = mtry_fc, 
             importance = "impurity",
             num.trees = 100,
             regularization.factor = .x, 
             regularization.usedepth = dep_par)

    })
}

boosted_svm <- mtry %>% 
  cross2(sim_new_sc) %>% 
  future_map(
    ~{
      map_boosted_svm(
        data = .x[[2]],
        mtry_fc = .x[[1]], 
        dep_par = 0
      )
    })

list(boosted_svm = boosted_svm) %>% saveRDS("rds/models/boosted_svm.Rds") 

#-------------------------------------------------------------------
# End of modeling
#-------------------------------------------------------------------