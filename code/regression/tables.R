# Package loading -------
library(tidyverse)

# Table 1. Variable Percentages -------------
# Mixture models 
df_mix <- read.table("auxiliar_results/df_mix.txt") 

print(
  df_mix %>% 
    mutate(lambdas_df = round(lambdas_df, 2)) %>% 
    filter(lambdas_df %in% c(0.1, 0.5, 0.9)) %>% 
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
                                 "Boosted SVM" = "boosted_svm",
                                 "GRRF" = "GRRF")) %>%

    xtable::xtable(), include.rownames = FALSE)

      

