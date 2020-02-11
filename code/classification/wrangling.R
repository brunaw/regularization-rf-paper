# Wrangling the results of the models 
library(tidyverse)

# Reading data ----------------
rds_names <- list.files("rds/classification", pattern = "_ws\\.rds$") 

rds <-  rds_names %>% 
  map(~{read_rds(paste0("rds/classification/", .x))})  

sizes <- map_dbl(
  .x = rds, 
  ~{
    res <- .x %>% map("result")
    if(is.null(res[[1]])){
      res <- .x
    }
    sum(!(res %>% map_lgl(is.null)))
  })

sizes # all ok!

# Manipulating ----------------
wrangle_rds <- map2(
  .x = rds, 
  .y = rds_names,
  ~{
    df <- data.frame(
      dataset = str_extract(.y, "[0-9]{1,2}") %>% as.numeric(),
      models = .x %>% map("model") %>% unlist(), 
      errors = .x %>% map("error") %>% unlist(),
      nvars = .x %>% map("nvars") %>% unlist()
    ) %>% 
      mutate(mtrys = rep(
        c(.x %>% map("mtrys") %>% unlist()), 4), 
        # Fixing the names
        models = as.character(models),
        # This is done because two type os models had their 
        # name switched when saving
        models = ifelse(models == "GRRF", "StdRF", 
                        ifelse(models == "Std_RF", "grrf", models))) 
    
    df
  }) %>% 
  bind_rows()

# Selecting optimal models and finding the mean number
# of variables and mean error in test set 

df  <- wrangle_rds %>% 
  mutate(sample = rep(rep(1:50, each = 20), 10)) %>% 
  group_by(dataset, models, sample) %>% 
  summarise(mean_error = mean(errors),
            sigma_error = sd(errors), 
            mean_vars = round(mean(nvars), 1)) %>% 
  group_by(dataset, models) %>% 
  arrange(dataset, models, mean_error, mean_vars) %>% 
  #filter(mean_error == min(mean_error)) %>% 
  slice(1) %>% 
  arrange(dataset, models, mean_vars, mean_error) %>% 
  slice(1) %>% 
  select(-sample) %>% 
  arrange(dataset, mean_error) %>% 
  mutate(mean_error = round(mean_error * 100, 2),
         sigma_error = round(sigma_error * 100, 2),
         full_error = 
           paste0(mean_error, " (", sigma_error, ")")) %>% 
  select(-mean_error, -sigma_error) %>% 
  gather(class, values, -dataset, -models) %>% 
  spread(models, values) %>% 
  ungroup() %>% 
  split(.$class) %>% 
  map(~{.x %>% 
      select(6, 5, 4, 3) %>% 
      setNames(c(
        "Standard RF", "GRRF", 
        "Boosted (RF)", "Boosted (MI)"))
  })

df %>% saveRDS("auxiliar_results/res_classification.rds")
#--------------------------------------------------------
