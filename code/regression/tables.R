# Package loading -------
library(tidyverse)

# Table 1. Variable Percentages -------------
# Mixture models 
df_mix <- read.table("auxiliar_results/df_mix.txt")

print(df_mix %>% 
        group_by(model, lambdas_df) %>% 
        summarise(mean_good = mean(good), 
                  mean_bad = mean(bad), 
                  mean_rmse = mean(rmse)) %>% 
        mutate_at(vars(mean_good, mean_bad), ~round(., 3)) %>% 
        mutate_at(vars(mean_good, mean_bad), 
                  ~scales::percent(., accuracy = 0.2 )) %>% 
        ungroup() %>% 
        mutate(model = 
                 forcats::fct_recode(model, 
                                     "Correlation" = "correlation", 
                                     "Boosted RF" = "boosted_rf", 
                                     "Boosted SVM" = "boosted_svm")) %>% 
        
        xtable::xtable(), include.rownames = FALSE)

#   $\text{Boosted}_{RF}$
# $\text{Boosted}_{SVM}$
      

# Table 2. Classification datasets info ------------
df_info <- read.table("auxiliar_results/dfs_info.txt") %>% 
        slice(-c(1, 8))  

print(df_info %>% 
              select(4, 1:3) %>% 
              mutate(dataset = str_replace(dataset, "\\.", " ")) %>% 
              setNames(c("Dataset", "Observations", "Features", "Classes")) %>% 
              xtable::xtable(), include.rownames = FALSE)

# Table 2. Classification datasets results ------------

res_classification <- readRDS("auxiliar_results/res_classification.rds")

print(res_classification$mean_vars %>% 
              mutate_all(as.numeric) %>% 
              mutate(
                      dataset = str_replace(df_info$dataset, "\\.", " "), 
                      columns = df_info$columns) %>% 
              mutate_at(vars(-columns, -dataset), funs(. /columns)) %>% 
              mutate_at(vars(-columns, -dataset), scales::percent) %>% 
              select(4, 1:3) %>% 
              bind_rows(
                      res_classification$mean_error %>% 
                              mutate(
                                      dataset = str_replace(df_info$dataset, "\\.", " "), 
                                      columns = df_info$columns) %>% 
                              select(4, 1:3)) %>% 
              xtable::xtable(), include.rownames = FALSE)

