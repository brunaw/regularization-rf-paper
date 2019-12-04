names <- list.files("data/classification/varselif", 
                    pattern = "data\\.txt$") %>% 
  .[-7] %>% 
  str_remove("\\.data\\.txt|\\.class\\.data\\.txt")

dfs_info <- final_dfs %>% 
  map_df(~{
    data.frame(
      rows = dim(.x)[1], 
      columns = dim(.x)[2], 
      classes = length(unique(.x$class)))
  }) %>% 
  mutate(dataset = names)

write.table(dfs_info, "data/classification/dfs_info.txt")