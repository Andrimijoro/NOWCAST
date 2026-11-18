
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

T <- as.Date("2026-03-30")
D <- 20L
H <- 26L
L <- 15L
B <- 1000L

# Minimum cohort age for delay-model estimation
M <- 52L

set.seed(123)

stopifnot(
  D >= 1,
  H >= 2,
  L <= D,
  M >= D
)

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
      (report_week - occurence_week) / 7
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

mature <- df %>%
  filter(
    occurence_week <= T - 7 * M
  )

stopifnot(nrow(mature) > 0)

# Pool event types with fewer than 30 mature events

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

# ------------------------------------------------------------
# Construct aggregated risk sets
#
# hits(d) = events reported exactly at d
# risk(d) = events with delay >= d
#
# For delay >20, observations contribute to the
# risk set through d=20 but never contribute a hit.
# ------------------------------------------------------------

risk_data <- mature %>%
  group_by(delay_group) %>%
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
    delay_group = factor(delay_group)
  )

# Verify all event types can be predicted

missing_groups <- setdiff(
  unique(type_lookup$delay_group),
  levels(risk_data$delay_group)
)

if (length(missing_groups) > 0) {
  stop(
    "No mature observations for: ",
    paste(missing_groups, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Complementary log-log GAM
#
# cloglog(pi_k(d)) = b_k(d)
# ------------------------------------------------------------

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
#plot(delay_model, pages = 1)

# ============================================================
# 4. REPORTING COMPLETENESS
# ============================================================

delay_grid <- tidyr::expand_grid(
  delay_group = levels(risk_data$delay_group),
  delay_week = 0:D
) %>%
  mutate(
    delay_group = factor(
      delay_group,
      levels = levels(risk_data$delay_group)
    )
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
  arrange(delay_group, delay_week) %>%
  group_by(delay_group) %>%
  mutate(
    rho = 1 - cumprod(1 - hazard)
  ) %>%
  ungroup() %>%
  mutate(
    delay_group = as.character(delay_group),
    rho = pmin(1, pmax(0, rho))
  )

# Recover original event types

type_cdf <- type_lookup %>%
  left_join(
    delay_cdf,
    by = "delay_group",
    relationship = "many-to-many"
  ) %>%
  select(
    event_type,
    delay_week,
    rho
  )

# Plot cumulative reporting distributions

ggplot(
  delay_cdf,
  aes(
    x = delay_week,
    y = rho,
    color = delay_group
  )
) +
  geom_line(linewidth = 1) +
  labs(
    title = "Cameroon",
    x = "Reporting delay (weeks)",
    y = "Cumulative reporting probability",
    color = "Event type"
  ) +
  theme_bw()

# ============================================================
# 5. HISTORICAL GAMMA-POISSON PRIOR
# ============================================================

# Historical: T-45,...,T-20
# Nowcasting: T,...,T-14

historical_weeks <- T - 7 * (D + 0:(H - 1L))

target_weeks <- T - 7 * (0:(L - 1L))

# Historical counts

historical_counts <- df %>%
  filter(
    occurence_week %in% historical_weeks
  ) %>%
  count(
    occurence_week,
    name = "N"
  )

# Add weeks with zero reported events

historical_panel <- tibble(
  occurence_week = historical_weeks
) %>%
  left_join(
    historical_counts,
    by = "occurence_week"
  ) %>%
  mutate(
    N = replace_na(N, 0L)
  )

# Method-of-moments empirical Bayes prior

mu <- mean(historical_panel$N)
v <- var(historical_panel$N)

stopifnot(mu > 0)

overdispersed <- v > mu

if (overdispersed) {
  
  alpha <- mu^2 / (v - mu)
  beta <- mu / (v - mu)
  
} else {
  
  # Fixed-intensity Poisson limiting model
  alpha <- NA_real_
  beta <- NA_real_
  
}

priors <- tibble(
  mu = mu,
  v = v,
  alpha = alpha,
  beta = beta,
  overdispersed = overdispersed
)

print(priors)

# ============================================================
# 6. HISTORICAL EVENT-TYPE WEIGHTS
# ============================================================

type_weights <- df %>%
  filter(
    occurence_week %in% historical_weeks
  ) %>%
  count(
    event_type,
    name = "n_type"
  ) %>%
  mutate(
    weight = n_type / sum(n_type)
  ) %>%
  select(event_type, weight)

# ============================================================
# 7. AGGREGATE REPORTING COMPLETENESS
# ============================================================

# rho_t(T) = sum_k w_k F_k(T-t)
#
# All quantities estimated using Cameroon only.

completeness <- tibble(
  occurence_week = target_weeks
) %>%
  mutate(
    age = as.integer(
      (T - occurence_week) / 7
    )
  ) %>%
  crossing(type_weights) %>%
  left_join(
    type_cdf,
    by = c(
      "event_type",
      "age" = "delay_week"
    )
  ) %>%
  group_by(
    occurence_week,
    age
  ) %>%
  summarise(
    rho = sum(weight * rho),
    .groups = "drop"
  )

stopifnot(
  all(is.finite(completeness$rho)),
  all(between(completeness$rho, 0, 1))
)

print(completeness)

# ============================================================
# 8. OBSERVED COUNTS AT VINTAGE T
# ============================================================

observed_counts <- df %>%
  filter(
    occurence_week %in% target_weeks
  ) %>%
  count(
    occurence_week,
    name = "Y"
  )

observed_panel <- tibble(
  occurence_week = target_weeks
) %>%
  left_join(
    observed_counts,
    by = "occurence_week"
  ) %>%
  mutate(
    Y = replace_na(Y, 0L)
  )

# ============================================================
# 9. GAMMA-POISSON NOWCAST
# ============================================================

# Prior:
# lambda_t ~ Gamma(alpha, beta)
#
# Likelihood:
# Y_t | lambda_t ~ Poisson(rho_t * lambda_t)
#
# Posterior:
# lambda_t | Y_t ~ Gamma(alpha + Y_t, beta + rho_t)

nowcast <- observed_panel %>%
  left_join(
    completeness,
    by = "occurence_week"
  ) %>%
  mutate(
    mu = mu,
    v = v,
    alpha = alpha,
    beta = beta,
    overdispersed = overdispersed,
    
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
    
    # Expected number of unreported events
    U_hat = if_else(
      overdispersed,
      (1 - rho) * alpha_post / beta_post,
      (1 - rho) * mu
    ),
    
    # Gamma-Poisson nowcast
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
  arrange(occurence_week)

# ============================================================
# 11. LATEST VINTAGE COUNTS
# ============================================================

latest_counts <- all_events %>%
  filter(
    occurence_week %in% target_weeks
  ) %>%
  count(
    occurence_week,
    name = "N_latest"
  )

nowcast_results <- nowcast %>%
  left_join(
    latest_counts,
    by = "occurence_week"
  ) %>%
  mutate(
    N_latest = replace_na(N_latest, 0L),
    country = country_name
  ) %>%
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
    upper,
    N_latest
  )

print(nowcast_results)

# ============================================================
# 12. PLOT
# ============================================================

plot_nowcast <- function() {
  
  ggplot(
    nowcast_results,
    aes(x = occurence_week)
  ) +
    
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      fill = "#356A87",
      alpha = 0.18
    ) +
    
    geom_line(
      aes(
        y = N_hat,
        color = "GP nowcast"
      ),
      linewidth = 0.9
    ) +
    
    geom_point(
      aes(
        y = N_hat,
        color = "GP nowcast"
      ),
      size = 2
    ) +
    
    geom_line(
      aes(
        y = Y,
        color = "At vintage T"
      ),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    
    geom_point(
      aes(
        y = Y,
        color = "At vintage T"
      ),
      size = 1.7
    ) +
    
    geom_line(
      aes(
        y = N_latest,
        color = "Latest vintage"
      ),
      linewidth = 0.9
    ) +
    
    geom_point(
      aes(
        y = N_latest,
        color = "Latest vintage"
      ),
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

plot_nowcast()

# ============================================================
# 13. EVALUATION
# ============================================================

evaluation <- nowcast_results %>%
  mutate(
    error_GP = N_hat - N_latest,
    error_offset = N_offset - N_latest,
    error_vintage = Y - N_latest
  )

evaluation %>%
  summarise(
    MSE_GP = mean(error_GP^2),
    MSE_offset = mean(error_offset^2, na.rm = TRUE),
    MSE_vintage = mean(error_vintage^2),
    coverage_95 = mean(
      N_latest >= lower &
        N_latest <= upper
    )
  )
