
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(mgcv)
library(lubridate)

# ============================================================
# 1. PARAMETERS
# ============================================================

T <- as.Date("2026-03-30")
D <- 20L
H <- 26L
L <- 15L
B <- 1000L

# Age required for fitting the delay distribution
M <- 52L

set.seed(123)

stopifnot(
  D >= 1,
  H >= 2,
  L >= 1,
  L <= D,
  M >= D
)

# ============================================================
# 2. PREPARE ACLED DATA
# ============================================================

# Weekly intervals are anchored at vintage T.
# Each week is (T-7, T], (T-14, T-7], etc.

week_T <- function(date) {
  T - 7 * floor(
    as.numeric(T - as.Date(date)) / 7
  )
}

# Prepare the full archive first.
# This is also used later to construct the latest-vintage
# reference counts.

all_events <- data %>%
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
  )

# Reconstruct the information available at vintage T

df <- all_events %>%
  filter(
    timestamp_date <= T,
    event_date <= T,
    delay_week >= 0
  )

stopifnot(
  nrow(df) > 0,
  all(
    as.numeric(df$report_week -
                 df$occurence_week) %% 7 == 0
  )
)

# IMPORTANT:
# Do not discard events with delay > D.
# They will be right-censored in the hazard model.

# ============================================================
# 3. ESTIMATE DISCRETE-TIME REPORTING HAZARDS
# ============================================================

# Use mature occurrence cohorts to reduce incomplete
# ascertainment of long reporting delays.

mature <- df %>%
  filter(
    occurence_week <= T - 7 * M
  )

stopifnot(nrow(mature) > 0)

# Group sparse event types

type_sizes <- mature %>%
  count(event_type, name = "n_type")

type_lookup <- df %>%
  distinct(event_type) %>%
  left_join(type_sizes, by = "event_type") %>%
  mutate(
    n_type = replace_na(n_type, 0L),
    delay_group = if_else(
      n_type >= 30,
      event_type,
      "Other"
    )
  ) %>%
  select(event_type, delay_group)

mature <- mature %>%
  left_join(type_lookup, by = "event_type")

# Check that every modeled group has mature observations

missing_groups <- setdiff(
  unique(type_lookup$delay_group),
  unique(mature$delay_group)
)

if (length(missing_groups) > 0) {
  stop(
    "No mature observations for delay group(s): ",
    paste(missing_groups, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Construct risk sets by country and event type
#
# hits(d) = number first reported at delay d
# risk(d) = number not reported before delay d
#
# Delays >20:
#   contribute to the risk set through d=20
#   but never contribute a hit.
#
# Exact delay =20:
#   contributes a hit at d=20.
# ------------------------------------------------------------

risk_data <- mature %>%
  group_by(country, delay_group) %>%
  group_modify(~ {
    
    x <- .x$delay_week
    
    tibble(
      delay_week = 0:D,
      
      hits = vapply(
        0:D,
        function(d) sum(x == d),
        integer(1)
      ),
      
      risk = vapply(
        0:D,
        function(d) sum(x >= d),
        integer(1)
      )
    )
  }) %>%
  ungroup() %>%
  filter(risk > 0) %>%
  mutate(
    country = factor(country),
    delay_group = factor(delay_group)
  )

# Summary of artificial right-censoring

censoring_summary <- mature %>%
  summarise(
    n_events = n(),
    n_censored = sum(delay_week > D),
    pct_censored = 100 * mean(delay_week > D)
  )

print(censoring_summary)

# ------------------------------------------------------------
# Discrete-time cloglog GAM
#
# cloglog(pi_{c,k}(d)) = b_k(d) + u_c
#
# u_c: country random intercept
# ------------------------------------------------------------

delay_model <- mgcv::gam(
  cbind(hits, risk - hits) ~
    delay_group +
    s(
      delay_week,
      by = delay_group,
      bs = "cr",
      k = 4
    ) +
    country,
  data = risk_data,
  family = binomial(link = "cloglog"),
  method = "REML"
)

summary(delay_model)
plot(delay_model, pages = 1)

# ============================================================
# 4. COUNTRY-SPECIFIC REPORTING COMPLETENESS
# ============================================================

# Predict hazard for every country, event type and delay.

delay_grid <- tidyr::expand_grid(
  country = levels(risk_data$country),
  delay_group = levels(risk_data$delay_group),
  delay_week = 0:D
) %>%
  mutate(
    country = factor(
      country,
      levels = levels(risk_data$country)
    ),
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

# Convert hazards to cumulative reporting probabilities.
#
# F_{c,k}(d) = 1 - prod_{j=0}^d (1 - pi_{c,k}(j))
#
# Do NOT force F(20)=1: delays >20 are right-censored.

delay_cdf <- delay_grid %>%
  arrange(
    country,
    delay_group,
    delay_week
  ) %>%
  group_by(country, delay_group) %>%
  mutate(
    rho = 1 - cumprod(1 - hazard)
  ) %>%
  ungroup() %>%
  mutate(
    country = as.character(country),
    delay_group = as.character(delay_group),
    rho = pmin(1, pmax(0, rho))
  ) %>%
  select(
    country,
    delay_group,
    delay_week,
    hazard,
    rho
  )

# Return to original event-type labels

type_cdf <- type_lookup %>%
  inner_join(
    delay_cdf,
    by = "delay_group",
    relationship = "many-to-many"
  )

# Visualize country-specific completeness curves

ggplot(
  delay_cdf,
  aes(x = delay_week, y = rho, group = country)
) +
  geom_line(alpha = 0.35) +
  facet_wrap(~delay_group) +
  labs(
    x = "Reporting delay (weeks)",
    y = "Cumulative reporting probability"
  ) +
  theme_bw()

# ============================================================
# 5. COUNTRY-SPECIFIC HISTORICAL PRIORS
# ============================================================

# Historical weeks T-45,...,T-20
# Target weeks T,T-1,...,T-14

historical_weeks <- T - 7 * (D + 0:(H - 1L))
target_weeks <- T - 7 * (0:(L - 1L))

# Countries represented in the fitted hazard model
countries <- levels(risk_data$country)

# Historical counts

historical_counts <- df %>%
  filter(
    occurence_week %in% historical_weeks,
    country %in% countries
  ) %>%
  count(
    country,
    occurence_week,
    name = "N"
  )

# Include weeks with zero observed events

historical_panel <- tidyr::expand_grid(
  country = countries,
  occurence_week = historical_weeks
) %>%
  left_join(
    historical_counts,
    by = c("country", "occurence_week")
  ) %>%
  mutate(
    N = replace_na(N, 0L)
  )

# Empirical-Bayes prior parameters
#
# lambda ~ Gamma(alpha, beta)
# N | lambda ~ Poisson(lambda)
#
# alpha = mu^2 / (v - mu)
# beta  = mu / (v - mu)

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

if (any(priors$mu <= 0)) {
  warning(
    "Countries with zero historical mean are omitted."
  )
}

priors <- priors %>%
  filter(mu > 0)

# ============================================================
# 6. HISTORICAL EVENT-TYPE WEIGHTS
# ============================================================

# w_{c,k}: historical proportion of event type k
# within country c.

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
  select(
    country,
    event_type,
    weight
  )

# ============================================================
# 7. COUNTRY-SPECIFIC AGGREGATE COMPLETENESS
# ============================================================

# rho_{c,t}(T) = sum_k w_{c,k} F_{c,k}(T-t)

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
      "country",
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
    occurence_week %in% target_weeks,
    country %in% priors$country
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

# Conditional on lambda:
#
# Y | lambda ~ Poisson(rho * lambda)
#
# Posterior:
#
# lambda | Y ~ Gamma(alpha + Y, beta + rho)
#
# E[U | Y] = (1-rho) * (alpha+Y)/(beta+rho)

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
    
    # Posterior parameters
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
    
    # Expected unreported events
    U_hat = if_else(
      overdispersed,
      (1 - rho) * alpha_post / beta_post,
      (1 - rho) * mu
    ),
    
    # GP nowcast
    N_hat = Y + U_hat,
    
    # Inverse-probability benchmark
    N_offset = if_else(
      rho > 0,
      Y / rho,
      NA_real_
    )
  )

# ============================================================
# 10. POSTERIOR PREDICTIVE INTERVALS
# ============================================================

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
    
    # Negative-binomial posterior predictive distribution
    U <- rnbinom(
      n = B,
      size = alpha_post,
      prob = beta_post /
        (beta_post + 1 - rho)
    )
    
  } else {
    
    # Limiting Poisson model
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
# 11. FINAL RESULTS AND LATEST VINTAGE
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

# Latest archive counts, including reports after T.
# Assumes 'data' contains the later archive.

latest_counts <- all_events %>%
  filter(
    occurence_week %in% target_weeks
  ) %>%
  count(
    country,
    occurence_week,
    name = "N_latest"
  )

nowcast_results <- nowcast_results %>%
  left_join(
    latest_counts,
    by = c("country", "occurence_week")
  ) %>%
  mutate(
    N_latest = replace_na(N_latest, 0L)
  )

print(nowcast_results)

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
    
    scale_x_date(
      date_labels = "%d %b"
    ) +
    
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

# Example
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
    aes(y = Y, color = "At vintage T"),
    linetype = "dashed"
  ) +
  geom_line(
    aes(y = N_latest, color = "Latest vintage")
  ) +
  facet_wrap(
    ~country,
    scales = "free_y"
  ) +
  scale_x_date(
    date_labels = "%d %b"
  ) +
  scale_color_manual(
    values = c(
      "GP nowcast" = "#235789",
      "At vintage T" = "#C45A22",
      "Latest vintage" = "#333333"
    )
  ) +
  labs(
    x = "Occurrence week",
    y = "Event count",
    color = NULL
  ) +
  theme_bw()
