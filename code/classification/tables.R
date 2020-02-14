# Package loading -------
library(tidyverse)

# Preparing Table 2
names <- list.files("data/classification/files/", 
                    pattern = "data\\.txt$") %>% 
        str_remove("\\.data\\.txt|\\.class\\.data\\.txt") %>% 
        str_remove("long\\.names") %>% 
        str_replace("\\.", " ") %>% 
        .[-7]

final_dfs <-  readRDS("data/classification/final_dfs.Rds")

dfs_info <- final_dfs %>% 
        map_df(~{
                data.frame(
                        rows = dim(.x)[1], 
                        columns = dim(.x)[2], 
                        classes = length(unique(.x$class)))
        }) %>% 
        mutate(dataset = names)

write.table(dfs_info, "data/classification/dfs_info.txt")

# Table 2. Classification datasets info ------------
df_info <- read.table("data/classification/dfs_info.txt") 
print(df_info %>% 
        select(4, 1:3) %>% 
        mutate(dataset = str_replace(dataset, "\\.", " ")) %>% 
        setNames(c("Dataset", "Observations", "Features", "Classes")) %>% 
        xtable::xtable(), include.rownames = FALSE)

# Table 3. Classification datasets results ------------
res_classification <- readRDS("auxiliar_results/res_classification.rds")
print(res_classification$mean_vars %>% 
        mutate_all(as.numeric) %>% 
        mutate_all(funs(round(./df_info$columns, 4)*100)) %>% 
        mutate(
          dataset = str_replace(df_info$dataset, "\\.", " ")) %>% 
        select(5, 1:4) %>% 
        bind_cols(
          res_classification$full_error) %>% 
        xtable::xtable(), include.rownames = FALSE)

