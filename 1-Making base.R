library(tidyverse)

####create the base data
##calculate the delay from this

df <- readRDS("data_df.rds")

data <- df %>%
  filter(event_date > as.Date("2024-06-30")) %>%
  group_by(id) %>%
  filter(date_id == min(date_id)) %>%
  mutate(interaction = as.character(interaction))

###add Daniel data
##get all of the folders where the event data are
folders <- list.dirs("data_daniel", recursive = FALSE) 


##create a function which reads the data inside each file and combines them.
read_one_folder <- function(folder_path) {
  files <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)
  
  folder_name <- basename(folder_path)
  
  files %>%
    lapply(function(f) {
      df <- readRDS(f) %>%
        filter(str_detect(region, "Africa")) %>%
        mutate(date_id = as.Date(folder_name))
      
      df
    }) %>%
    bind_rows()
}

##combine the june-august data

data_Daniel <- folders %>%
  lapply(read_one_folder) %>%
  bind_rows()

####clean Daniel data


df <- data_Daniel %>%
  rename(id = event_id_cnty) %>%
  mutate(timestamp_date = as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC")) %>%
  select(colnames(data)) %>%
  mutate(event_date = as.Date(event_date), interaction = as.character(interaction)) %>%
  filter(event_date > as.Date("2024-06-30"), !(id %in% unique(data$id))) %>%
  group_by(id) %>%
  filter(date_id == min(date_id))

####BASE DATA

df$latitude <- as.numeric(df$latitude)
df$longitude <- as.numeric(df$longitude)

base_df <- rbind(data, df)

saveRDS(base_df, "base_df.rds")

