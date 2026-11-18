####require the data from 3- Analysis Data


library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(mgcv)
library(lubridate)

# ============================================================
# 1. PARAMETERS
# ============================================================

T <- as.Date("2026-03-23")+7  # Historical vintage
D <- 20L                   # Maximum reporting delay (weeks)
H <- 26L                   # Historical prior window
L <- 15L                   # Number of weeks to nowcast
B <- 1000L                 # Posterior predictive simulations

set.seed(123)

stopifnot(
  D == 20L,
  H == 26L,
  L <= D,
  lubridate::wday(T, week_start = 1) == 1L
)

# ============================================================
# 2. PREPARE ACLED DATA
# ============================================================


# Define weekly intervals relative to T
week_T <- function(date) {
  T - 7 * floor(as.numeric(T - as.Date(date)) / 7)
}

data <- readRDS("base_df.rds")

df <- data %>%
  mutate(
    event_date = as.Date(event_date),
    timestamp_date = as.Date(timestamp_date),
    event_type = as.character(event_type),
    country = as.character(country)
  ) %>%
  filter(
    !is.na(id),
    !is.na(country),
    !is.na(event_type),
    !is.na(event_date),
    !is.na(timestamp_date),
    timestamp_date <= T,
    event_date <= T,
    timestamp_date >= event_date
  ) %>%
  arrange(timestamp_date) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    occurence_week = week_T(event_date),
    report_week = week_T(timestamp_date),
    delay_week = as.integer(
      (report_week - occurence_week) / 7
    )
  ) %>%
  filter(
    !is.na(delay_week),
    delay_week >= 0
  )

# Check whether dates genuinely represent full weeks
stopifnot(
  all(
    as.numeric(df$report_week - df$occurence_week) %% 7 == 0
  )
)

# Apply the manuscript's D-week reporting horizon.
# Events with longer delays are excluded from estimation,
# not silently reclassified as having delay D.
if (any(df$delay_week > D)) {
  warning(
    "Some delays exceed D=20. These events are excluded. ",
    "Check the bounded-delay assumption."
  )
}

df <- df %>%
  filter(delay_week <= D)

# ============================================================
# 3. ESTIMATE DISCRETE-TIME REPORTING HAZARDS
# ============================================================

# Use mature occurrence cohorts to estimate delay hazards.
mature <- df %>%
  filter(occurence_week <= T - 7 * D)

# Pool sparse event types for hazard estimation.
# Main event types retain their own delay curves.
type_sizes <- mature %>%
  count(event_type, name = "n_type")

type_lookup <- type_sizes %>%
  mutate(
    delay_group = if_else(
      n_type >= 30,
      event_type,
      "Other"
    )
  ) %>%
  select(event_type, delay_group)

mature <- mature %>%
  left_join(type_lookup, by = "event_type")

# Risk-set construction:
# risk(d) = number of events with delay >= d
# hits(d) = number of events with delay = d

risk_data <- mature %>%
  group_by(delay_group) %>%
  group_modify(~ {
    x <- .x$delay_week
    
    tibble(
      delay_week = 0:(D - 1L),
      hits = vapply(
        0:(D - 1L),
        function(d) sum(x == d),
        integer(1)
      ),
      risk = vapply(
        0:(D - 1L),
        function(d) sum(x >= d),
        integer(1)
      )
    )
  }) %>%
  ungroup() %>%
  filter(risk > 0) %>%
  mutate(
    delay_group = factor(delay_group)
  )

# Complementary log-log GAM
#
# cloglog(hazard_k(d)) = b_k(d)

delay_model <- mgcv::gam(
  cbind(hits, risk - hits) ~
    delay_group +
    s(
      delay_week,
      by = delay_group,
      bs = "cr",
      k = 4
    ),
  data = risk_data,
  family = binomial(link = "cloglog"),
  method = "REML"
)

summary(delay_model)

# ============================================================
# 4. REPORTING COMPLETENESS F_k(d)
# ============================================================

delay_grid <- tidyr::expand_grid(
  delay_group = levels(risk_data$delay_group),
  delay_week = 0:(D - 1L)
) %>%
  mutate(
    delay_group = factor(
      delay_group,
      levels = levels(risk_data$delay_group)
    )
  )

delay_grid$hazard <- predict(
  delay_model,
  newdata = delay_grid,
  type = "response"
)

# Convert hazards into cumulative reporting probabilities
delay_cdf <- delay_grid %>%
  arrange(delay_group, delay_week) %>%
  group_by(delay_group) %>%
  mutate(
    rho = 1 - cumprod(1 - hazard)
  ) %>%
  ungroup() %>%
  select(delay_group, delay_week, rho) %>%
  bind_rows(
    tibble(
      delay_group = factor(
        levels(risk_data$delay_group),
        levels = levels(risk_data$delay_group)
      ),
      delay_week = D,
      rho = 1
    )
  ) %>%
  mutate(
    delay_group = as.character(delay_group),
    rho = pmin(1, pmax(0, rho))
  )

# Completion probability by actual event type
type_cdf <- type_lookup %>%
  left_join(
    delay_cdf,
    by = "delay_group",
    relationship = "many-to-many"
  )

# ============================================================
# 5. COUNTRY-SPECIFIC HISTORICAL PRIORS
# ============================================================

# Historical weeks: T-45,...,T-20 (26 weeks)

historical_weeks <- T - 7 * (20:45)
target_weeks <- max(df$occurence_week) - 7 * (0:(L - 1L))

countries <- sort(unique(df$country))
event_types <- sort(unique(df$event_type))

# Explicit zero counts are essential.
# Otherwise historical means will be biased upward.

historical_counts <- df %>%
  filter(occurence_week %in% historical_weeks) %>%
  count(
    country,
    occurence_week,
    name = "N"
  )

historical_panel <- tidyr::expand_grid(
  country = countries,
  occurence_week = historical_weeks
) %>%
  left_join(
    historical_counts,
    by = c("country", "occurence_week")
  ) %>%
  mutate(N = replace_na(N, 0L))

# Historical mean, variance, Gamma prior parameters
priors <- historical_panel %>%
  group_by(country) %>%
  summarise(
    mu = mean(N),
    v = var(N),
    .groups = "drop"
  ) %>%
  mutate(
    overdispersed = v > mu,
    alpha = if_else(
      overdispersed,
      mu^2 / (v - mu),
      NA_real_
    ),
    beta = if_else(
      overdispersed,
      mu / (v - mu),
      NA_real_
    )
  )

# The paper assumes mu > 0.
# Do not silently manufacture a positive prior.
if (any(priors$mu <= 0)) {
  warning(
    "Countries with zero historical mean will be omitted."
  )
}

priors <- priors %>%
  filter(mu > 0)

# ============================================================
# 6. HISTORICAL EVENT-TYPE WEIGHTS
# ============================================================

# Estimate w_{c,k} from the same 26-week historical window.

type_weights <- df %>%
  filter(
    occurence_week %in% historical_weeks,
    country %in% priors$country
  ) %>%
  count(
    country,
    event_type,
    name = "n_type"
  ) %>%
  group_by(country) %>%
  mutate(
    weight = n_type / sum(n_type)
  ) %>%
  ungroup() %>%
  select(country, event_type, weight)

# ============================================================
# 7. CALCULATE REPORTING COMPLETENESS
# ============================================================

# rho_{c,t}(T) = sum_k w_{c,k} F_k(T-t)

completeness <- tidyr::expand_grid(
  country = priors$country,
  occurence_week = target_weeks
) %>%
  mutate(
    age = as.integer(
      (T - occurence_week) / 7
    )
  ) %>%
  left_join(
    type_weights,
    by = "country",
    relationship = "many-to-many"
  ) %>%
  left_join(
    type_cdf,
    by = c(
      "event_type",
      "age" = "delay_week"
    )
  ) %>%
  group_by(
    country,
    occurence_week,
    age
  ) %>%
  summarise(
    rho = sum(weight * rho),
    .groups = "drop"
  )

stopifnot(
  all(is.finite(completeness$rho)),
  all(completeness$rho >= 0),
  all(completeness$rho <= 1)
)

# ============================================================
# 8. OBSERVED COUNTS AT VINTAGE T
# ============================================================

observed_counts <- df %>%
  filter(
    occurence_week %in% target_weeks
  ) %>%
  count(
    country,
    occurence_week,
    name = "Y"
  )

observed_panel <- tidyr::expand_grid(
  country = priors$country,
  occurence_week = target_weeks
) %>%
  left_join(
    observed_counts,
    by = c("country", "occurence_week")
  ) %>%
  mutate(
    Y = replace_na(Y, 0L)
  )

# ============================================================
# 9. GAMMA-POISSON POSTERIOR NOWCAST
# ============================================================

nowcast <- observed_panel %>%
  left_join(
    completeness,
    by = c("country", "occurence_week")
  ) %>%
  left_join(
    priors,
    by = "country"
  ) %>%
  mutate(
    # Posterior Gamma parameters
    alpha_post = if_else(
      overdispersed,
      alpha + Y,
      NA_real_
    ),
    beta_post = if_else(
      overdispersed,
      beta + rho,
      NA_real_
    ),
    
    # Posterior expectation of unreported events
    U_hat = if_else(
      overdispersed,
      (1 - rho) * alpha_post / beta_post,
      (1 - rho) * mu
    ),
    
    # Posterior mean of total event counts
    N_hat = Y + U_hat,
    
    # Simple inverse-probability benchmark
    N_offset = if_else(
      rho > 0,
      Y / rho,
      NA_real_
    )
  )

# ============================================================
# 10. POSTERIOR PREDICTIVE INTERVALS
# ============================================================

# Sample U directly from the posterior predictive NB
# distribution. If the historical variance <= mean,
# use the limiting fixed-intensity Poisson model.

simulate_total <- function(
    Y, rho, mu,
    overdispersed,
    alpha_post, beta_post,
    B = 1000L
) {
  
  if (rho >= 1) {
    return(rep(Y, B))
  }
  
  if (overdispersed) {
    
    U <- rnbinom(
      n = B,
      size = alpha_post,
      prob = beta_post /
        (beta_post + 1 - rho)
    )
    
  } else {
    
    U <- rpois(
      n = B,
      lambda = (1 - rho) * mu
    )
  }
  
  Y + U
}

nowcast <- nowcast %>%
  mutate(
    draws = pmap(
      list(
        Y,
        rho,
        mu,
        overdispersed,
        alpha_post,
        beta_post
      ),
      ~ simulate_total(
        Y = ..1,
        rho = ..2,
        mu = ..3,
        overdispersed = ..4,
        alpha_post = ..5,
        beta_post = ..6,
        B = B
      )
    ),
    lower = map_dbl(
      draws,
      ~ unname(quantile(.x, 0.025))
    ),
    upper = map_dbl(
      draws,
      ~ unname(quantile(.x, 0.975))
    )
  ) %>%
  arrange(country, occurence_week)

# ============================================================
# 11. FINAL RESULTS
# ============================================================

nowcast_results <- nowcast %>%
  select(
    country,
    occurence_week,
    age,
    Y,
    rho,
    mu,
    v,
    alpha,
    beta,
    U_hat,
    N_hat,
    N_offset,
    lower,
    upper
  )

print(nowcast_results)

# Counts at the latest vintage
latest_counts <- data %>%
  mutate(
    occurence_week = week_T(event_date)
  ) %>%
  filter(occurence_week %in% target_weeks) %>%
  distinct(id, .keep_all = TRUE) %>%
  count(country, occurence_week, name = "N_latest")

nowcast_results <- nowcast_results %>%
  left_join(
    latest_counts,
    by = c("country", "occurence_week")
  ) %>%
  mutate(N_latest = replace_na(N_latest, 0L))


# ============================================================
# 12. VISUALIZATION
# ============================================================

plot_nowcast <- function(country_name) {
  
  pdat <- nowcast_results %>%
    filter(country == country_name) %>%
    arrange(occurence_week)
  
  ggplot(pdat, aes(x = occurence_week)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "#356A87",
      alpha = 0.18
    ) +
    geom_line(
      aes(y = N_hat, color = "GP nowcast"),
      linewidth = 0.9
    ) +
    geom_point(
      aes(y = N_hat, color = "GP nowcast"),
      size = 2
    ) +
    geom_line(
      aes(y = Y, color = "At vintage T"),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    geom_point(
      aes(y = Y, color = "At vintage T"),
      size = 1.7
    ) +
    geom_line(
      aes(y = N_latest, color = "Latest vintage"),
      linewidth = 0.9
    ) +
    geom_point(
      aes(y = N_latest, color = "Latest vintage"),
      size = 1.7
    ) +
    scale_color_manual(
      values = c(
        "GP nowcast" = "#235789",
        "At vintage T" = "#C45A22",
        "Latest vintage" = "#333333"
      )
    ) +
    scale_x_date(date_labels = "%d %b") +
    labs(
      title = country_name,
      x = "Occurrence week",
      y = "Number of conflict events",
      color = NULL,
      caption = "Shaded area: 95% posterior predictive interval"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

# Example:
plot_nowcast("Cameroon")

# All countries
ggplot(
  nowcast_results,
  aes(x = occurence_week)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#356A87",
    alpha = 0.15
  ) +
  geom_line(
    aes(y = N_hat, color = "GP nowcast")
  ) +
  geom_line(
    aes(y = Y, color = "Observed"),
    linetype = "dashed"
  ) +
  facet_wrap(~country, scales = "free_y") +
  scale_x_date(date_labels = "%d %b") +
  labs(
    x = "Occurrence week",
    y = "Event count",
    color = NULL
  ) +
  theme_bw()
