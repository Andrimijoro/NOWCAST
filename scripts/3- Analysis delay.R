library(tidyverse)
library(slider)
library(scales)

data <- readRDS("base_df.rds")

data <- data %>%
  mutate(timestamp_date = as.Date(as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC")),
         delay = floor(as.numeric(timestamp_date-event_date)/7),
         occurence_week = floor_date(event_date, unit = "week", week_start = 1),
         report_week = floor_date(timestamp_date, unit = "week", week_start = 1)) %>%
  mutate(delay = ifelse(delay>20,20,delay))

###plot delay median

#####another approach 

df1 <- data %>%
  group_by(occurence_week) %>%
  summarise(
    n = n(),
    median = median(delay, na.rm = TRUE),
    q25 = quantile(delay, 0.25, na.rm = TRUE),
    q75 = quantile(delay, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(occurence_week <= max(occurence_week) - 20*7)

p <- df1 %>%
  pivot_longer(
    cols = c(median, q25, q75),
    names_to = "Statistic",
    values_to = "Delay"
  ) %>%
  mutate(
    Statistic = factor(
      Statistic,
      levels = c("median", "q75", "q25"),
      labels = c("Median", "75th percentile", "25th percentile")
    )
  ) %>%
  ggplot(aes(x = occurence_week, y = Delay,
             color = Statistic, linetype = Statistic)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    values = c(
      "Median"          = "#0072B2",  # Deep blue
      "75th percentile" = "#D55E00",  # Vermilion
      "25th percentile" = "#009E73"   # Teal green
    )
  ) +
  scale_linetype_manual(
    values = c("Median" = "solid",
               "75th percentile" = "dashed",
               "25th percentile" = "dashed")
  ) +
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b %Y"
  ) +
  labs(
    x = "Occurrence date",
    y = "Reporting delay (weeks)",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_line(
      color = "grey90", linewidth = 0.3
    ),
    axis.text = element_text(color = "grey20"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = grid::unit(1.5, "cm")
  )

p


####plot with report date

df2 <- data %>%
  group_by(report_week) %>%
  mutate(median = median(delay), 
         q25 = quantile(delay, 0.25),
         q75 = quantile(delay, 0.75)) %>%
  distinct(report_week, median, q25, q75) %>%
  ungroup() %>%
  filter(report_week >= min(report_week) + 20*7)

p1 <- df2 %>%
  pivot_longer(
    cols = c(median, q25, q75),
    names_to = "Statistic",
    values_to = "Delay"
  ) %>%
  mutate(
    Statistic = factor(
      Statistic,
      levels = c("median", "q25", "q75"),
      labels = c("Median", "75th percentile", "25th percentile")
    )
  ) %>%
  ggplot(aes(x = report_week, y = Delay,
             color = Statistic, linetype = Statistic)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    values = c(
      "Median"          = "#0072B2",  # Deep blue
      "75th percentile" = "#D55E00",  # Vermilion
      "25th percentile" = "#009E73"   # Teal green
    )
  ) +
  scale_linetype_manual(
    values = c("Median" = "solid",
               "75th percentile" = "dashed",
               "25th percentile" = "dashed")
  ) +
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b %Y"
  ) +
  labs(
    x = "Report date",
    y = "Reporting delay (weeks)",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_line(
      color = "grey90", linewidth = 0.3
    ),
    axis.text = element_text(color = "grey20"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = grid::unit(1.5, "cm")
  )

p1








library(dplyr)
library(tidyr)
library(ggplot2)

df1 <- data %>%
  group_by(occurence_week) %>%
  summarise(
    n = n(),
    median = median(delay, na.rm = TRUE),
    q25 = quantile(delay, 0.25, na.rm = TRUE),
    q75 = quantile(delay, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    occurence_week <= max(occurence_week) - 20 * 7
  )

# Scaling factor for the secondary axis
scale_factor <- max(df1$n, na.rm = TRUE) /
  max(df1$q75, df1$median, df1$q25, na.rm = TRUE)

df_long <- df1 %>%
  pivot_longer(
    cols = c(median, q25, q75),
    names_to = "Statistic",
    values_to = "Delay"
  ) %>%
  mutate(
    Statistic = factor(
      Statistic,
      levels = c("median", "q75", "q25"),
      labels = c("Median", "75th percentile", "25th percentile")
    )
  )

p <- ggplot() +
  
  # Unsmoothed delay quantiles
  geom_line(
    data = df_long,
    aes(
      x = occurence_week,
      y = Delay,
      color = Statistic,
      linetype = Statistic
    ),
    linewidth = 0.9
  ) +
  
  # Smoothed weekly event counts only
  geom_smooth(
    data = df1,
    aes(
      x = occurence_week,
      y = n / scale_factor,
      color = "Number of events",
      linetype = "Number of events"
    ),
    method = "loess",
    formula = y ~ x,
    span = 0.1,
    se = FALSE,
    linewidth = 1
  ) +
  
  scale_color_manual(
    values = c(
      "Median" = "#0072B2",
      "75th percentile" = "#D55E00",
      "25th percentile" = "#009E73",
      "Number of events" = "#555555"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "Median" = "solid",
      "75th percentile" = "dashed",
      "25th percentile" = "dashed",
      "Number of events" = "solid"
    )
  ) +
  
  scale_y_continuous(
    name = "Reporting delay (weeks)",
    sec.axis = sec_axis(
      ~ . * scale_factor,
      name = "Number of events"
    )
  ) +
  
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b %Y"
  ) +
  
  labs(
    x = "Occurrence date",
    color = NULL,
    linetype = NULL
  ) +
  
  theme_classic(
    base_size = 11,
    base_family = "serif"
  ) +
  
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    axis.text = element_text(color = "grey20"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.title.y.right = element_text(color = "#555555"),
    legend.key.width = grid::unit(1.5, "cm")
  )

p

ggsave(
  "reporting_delay_and_count.pdf",
  plot = p,
  width = 11,
  height = 6,
  device = "pdf",
  bg = "white"
)