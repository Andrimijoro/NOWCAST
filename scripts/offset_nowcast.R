#####This is a script for the offset model

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(mgcv)
library(lubridate)

# ============================================================
# 1. PARAMETERS
# ============================================================

country_name <- "Sudan"

T <- as.Date("2025-03-29")
D <- 21L
L <- 20L
B <- 1000L

# Minimum cohort age for delay-model estimation
M <- 25L

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

# Cameroon only
# Keep the full archive for later-vintage evaluation

data <- readRDS("base_df.rds")

all_events <- data %>%
  filter(country == country_name) %>%
  mutate(
    event_date = as.Date(event_date),
    timestamp_date = as.Date(timestamp_date),
    event_type = as.character(event_type),
    country = as.character(country)
  ) %>%
  filter(
    !is.na(id),
    !is.na(event_date),
    !is.na(timestamp_date),
    !is.na(event_type),
    timestamp_date >= event_date
  ) %>%
  arrange(timestamp_date) %>%
  distinct(id, .keep_all = TRUE) %>%
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
    date_id <= T,
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

mature <- df %>%
  filter(
    occurence_week <= T - 7 * M
  )

stopifnot(nrow(mature) > 0)

x <- mature$delay_week
x <- x[!is.na(x)]

risk_data <- tibble(
  delay_week = 0:(D - 1L),
  hits = vapply(0:(D - 1L), function(d) sum(x == d), integer(1)),
  risk = vapply(0:(D - 1L), function(d) sum(x >= d), integer(1))
) %>%
  filter(risk > 0)

# Pool event types with fewer than 30 mature events

# ------------------------------------------------------------
# Complementary log-log GAM
#
# cloglog(pi_k(d)) = b_k(d)
# ------------------------------------------------------------

delay_model <- mgcv::gam(
  cbind(hits, risk - hits) ~
    s(
      delay_week,
      bs = "tp"
    ),
  data = risk_data,
  family = binomial(link = "cloglog"),
  method = "REML"
)

summary(delay_model)
plot(delay_model, pages = 1)

# ============================================================
# 4. REPORTING COMPLETENESS
# ============================================================

delay_grid <- tidyr::expand_grid(
  delay_week = 0:D
)

# Predicted reporting hazards

delay_grid$hazard <- predict(
  delay_model,
  newdata = delay_grid,
  type = "response"
)

# Cumulative reporting probabilities:
#
# F_k(d) = 1 - prod_{j=0}^d (1 - pi_k(j))

delay_cdf <- delay_grid %>%
  arrange(delay_week) %>%
  mutate(
    rho = 1 - cumprod(1 - hazard)
  ) 

