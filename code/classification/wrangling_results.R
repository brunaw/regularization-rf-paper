# Code for wrangling the results of the models rerun 10 times
# Bruna Wundervald, November 2019 
library(tidyverse)

# Reading data ----------------
rds_names <- list.files("rds", pattern = "10\\.rds$") %>% 
  .[-c(1, 9)]

# Datasets 1 and 8 still have issues (?)

rds <-  rds_names %>% 
  map(~{read_rds(paste0("rds/", .x))})  

sizes <- map_dbl(
  .x = rds, 
  ~{
    res <- .x %>% map("result")
    if(is.null(res[[1]])){
      res <- .x
    }
    sum(!(res %>% map_lgl(is.null)))
  })



# Manipulating ----------------
wrangle_rds10 <- map2(
  .x = rds, 
  .y = rds_names,
  ~{
    res <- .x %>% map("result")
    #if(is.null(res[[1]])){
    #  res <- .x
    #}
    correct <- sum(!(res %>% map_lgl(is.null)))
    
    models <- res %>% map("model") %>% unlist()
    errors <- res %>% map("error") %>% unlist()
    nvars <- res %>% map("nvars") %>% unlist()
    mtrys <- c(NA, 
               rep(res %>% map("mtrys") %>% unlist() %>% 
                     unique(), 4))
    
    df <- data.frame(
      dataset = str_extract(.y, "[0-9]{1}"),
      models = models, 
      errors = errors,
      nvars = nvars
    ) %>% 
      mutate(mtrys = rep(mtrys, correct))
    
    df
  }) %>% 
  bind_rows()

# Selecting optimal models and finding the mean number
# of variables and mean error in test set 

df  <- wrangle_rds10 %>% 
  filter(models %in% c("MI", "GRRF", "RF")) %>% 
  group_by(dataset, models, mtrys) %>% 
  summarise(mean_error = median(errors),  # Median, not mean
            mean_vars = round(median(nvars), 1)) %>% 
  group_by(dataset, models) %>% 
  filter(mean_error == min(mean_error)) %>% 
  mutate(mean_error = scales::percent(mean_error)) %>% 
  arrange(mean_vars) %>% 
  slice(1) %>% 
  select(-mtrys) %>% 
  ungroup() %>% 
  gather(class, values, -dataset, -models) %>% 
  spread(models, values) %>% 
  split(.$class) %>% 
  map(~{.x %>% 
      select(5, 4, 3) %>% 
      setNames(c("Standard RF", 
                 "Boosted (RF)", "Boosted (MI)"))
  }) 

df %>% saveRDS("auxiliar_results/res_classification.rds")

