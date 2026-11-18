####This updates the base data 
###requires base_df.rds

library(tidyverse)

base_df <- readRDS("data/base_df.rds")
date_update <- max(base_df$date_id)


for (i in 1:60) {

  date_update <- date_update+7
  df <- readRDS(paste0("data/data_", date_update, ".rds"))
  
  df <- df %>%
    rename(id = event_id_cnty) %>%
    select(colnames(base_df)) %>%
    mutate(event_date = as.Date(event_date), interaction = as.character(interaction),
           longitude = as.numeric(longitude), latitude = as.numeric(latitude)) %>%
    filter(event_date > as.Date("2024-06-30"), !(id %in% unique(base_df$id))) 
  
  base_df <- rbind(base_df, df)
  print(paste(i, date_update, nrow(df)))
  
}


saveRDS(base_df, "base_df.rds")

