#####This is a script for the offset model

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(mgcv)
library(lubridate)
library(MASS)

# ============================================================
# 1. PARAMETERS
# ============================================================

country_name <- "Somalia"
data <- readRDS("base_df.rds")
data <- filter(data, country == country_name)

#this is the last mature t
T_end <- max(data$date_id) - 20*7 + 2

#this is the first analysis date, a Tuesday.
T_start <- T_end - 9*7

# Minimum cohort age for delay-model estimation
M <- 20L
D <- 20L
W <- 26
###number of window
w <- 10

pdf("all_plots.pdf", width = 10, height = 6)

for (i in 1:w){
  
  print(i)
  T <- T_start + (i-1)*7###the download date
  set.seed(123)
  
  # ============================================================
  # 2. PREPARE DATA
  # ============================================================
  
  # Define weekly intervals relative to T
  week_T <- function(date) {
    T - 7 * floor(
      as.numeric(T - as.Date(date)) / 7
    )
  }
  
  
  all_events <- data %>%
    mutate(
      event_date = as.Date(event_date),
      timestamp_date = as.Date(timestamp_date),
      event_type = as.character(event_type),
      country = as.character(country)
    ) %>%
    arrange(timestamp_date) %>%
    mutate(
      occurence_week = week_T(event_date),
      report_week = week_T(timestamp_date),
      delay_week = as.integer(
        (timestamp_date - event_date) / 7
      )
    )
  
  # Reconstruct vintage T
  
  df <- all_events %>%
    filter(
      timestamp_date <= T,
      event_date <= T,
      delay_week >= 0
    )
  
  stopifnot(nrow(df) > 0)
  
  # Do not remove delays >20.
  # They contribute right-censored observations.
  
  # ============================================================
  # 3. REPORTING DELAY DISTRIBUTION
  # ============================================================
  
  # Mature occurrence cohorts
  
  # ============================================================
  # 1. CONSTRUCT REPORTING TRIANGLE
  # ============================================================
  
  regression_table <- df %>%
    filter(
      occurence_week >= T - 14 * D,
      occurence_week <= T,
      between(delay_week, 0L, D)
    ) %>%
    count(occurence_week, delay_week, name = "n")
  
  reporting_triangle <- expand_grid(
    occurence_week = seq(T - 14 * D, T, by = "week"),
    delay_week = 0:D
  ) %>%
    mutate(
      age = as.integer(T - occurence_week) %/% 7,
      t = dense_rank(occurence_week),
      observed = delay_week <= age
    ) %>%
    left_join(
      regression_table,
      by = c("occurence_week", "delay_week")
    ) %>%
    mutate(
      n = if_else(
        observed,
        coalesce(n, 0L),
        NA_integer_
      )
    )
  
  observed_data <- reporting_triangle %>%
    filter(observed)
  
  prediction_data <- reporting_triangle %>%
    filter(!observed)
  
  
  # ============================================================
  # 2. ESTIMATE JOINT MODEL
  # ============================================================
  
  joint_model <- mgcv::gam(
    n ~ s(t, k = 4) +
      s(t, delay_week) ,
    family = nb(link = "log"),
    method = "REML",
    data = observed_data
  )
  
  
  # ============================================================
  # 3. PREDICT UNREPORTED COUNTS
  # ============================================================
  
  prediction_data <- prediction_data %>%
    mutate(
      predicted = predict(
        joint_model,
        newdata = prediction_data,
        type = "response"
      )
    )
  
  
  # ============================================================
  # 4. COMPUTE NOWCAST
  # ============================================================
  
  reported_data <- observed_data %>%
    group_by(occurence_week) %>%
    summarise(
      reported = sum(n),
      .groups = "drop"
    )
  
  unreported_data <- prediction_data %>%
    group_by(occurence_week) %>%
    summarise(
      unreported = sum(predicted),
      .groups = "drop"
    )
  
  nowcast_data <- reported_data %>%
    left_join(
      unreported_data,
      by = "occurence_week"
    ) %>%
    mutate(
      unreported = coalesce(unreported, 0),
      nowcast = reported + unreported
    ) %>%
    arrange(occurence_week)
  
  
  # ============================================================
  # 5. EVALUATION DATA
  # ============================================================
  
  true_counts <- all_events %>%
    filter(
      occurence_week >= T - 7 * D,
      occurence_week <= T_end,
      between(delay_week, 0L, D)
    ) %>%
    count(occurence_week, name = "true_n")
  
  data_final <- nowcast_data %>%
    filter(
      occurence_week >= T - 7 * D,
      occurence_week <= T_end
    ) %>%
    left_join(
      true_counts,
      by = "occurence_week"
    ) %>%
    mutate(
      true_n = coalesce(true_n, 0L),
      error = nowcast - true_n
    )
  
  p <- ggplot(data_final, aes(x = occurence_week)) +
    geom_line(aes(y = true_n, colour = "Final count")) +
    geom_line(aes(y = nowcast, colour = "Nowcast")) +
    geom_point(aes(y = reported, colour = "Reported")) +
    scale_colour_manual(
      values = c(
        "Final count" = "black",
        "Nowcast" = "red",
        "Reported" = "blue"
      )
    ) +
    labs(
      x = "Occurrence week",
      y = "Event count",
      colour = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  
}

dev.off()
