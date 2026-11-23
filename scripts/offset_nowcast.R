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

pdf("all_plots.pdf", width = 10, height = 6)

for (i in 1:10){
  
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
  
  mature <- df %>%
    filter(occurence_week <= T-M*7 & occurence_week >= T-(D+W-1)*7) 
  
  stopifnot(nrow(mature) > 0)
  
  x <- mature$delay_week
  x <- x[!is.na(x)]
  
  risk_data <- tibble(
    delay_week = 0:D,
    hits = vapply(0:D, function(d) sum(x == d), integer(1)),
    risk = vapply(0:D, function(d) sum(x >= d), integer(1))
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
  
  #summary(delay_model)
  #plot(delay_model, pages = 1)
  
  # ============================================================
  # 4. REPORTING COMPLETENESS
  # ============================================================
  
  delay_grid <- tidyr::expand_grid(
    delay_week = 0:D
  )
  
  delay_grid <- delay_grid %>%
    arrange(delay_week)
  
  # Model matrix and fitted hazards
  X <- predict(delay_model, newdata = delay_grid, type = "lpmatrix")
  eta <- drop(X %*% coef(delay_model))
  hazard <- delay_model$family$linkinv(eta)
  
  # Cumulative reporting probability
  S <- cumprod(1 - hazard)
  rho <- 1 - S
  
  # Delta-method gradient
  w <- delay_model$family$mu.eta(eta) / (1 - hazard)
  G <- apply(sweep(X, 1, w, "*"), 2, cumsum)
  G <- sweep(G, 1, S, "*")
  
  # Asymptotic standard error
  V <- vcov(delay_model, unconditional = TRUE)
  se <- sqrt(rowSums((G %*% V) * G))
  
  # 95% confidence intervals
  delay_cdf <- delay_grid %>%
    mutate(
      hazard = hazard,
      rho = rho,
      rho_lower = pmax(0, rho - 2 * se),
      rho_upper = pmin(1, rho + 2 * se)
    )
  
  ###Table 
  
  regression_table <- df %>%
    filter(occurence_week >= T - D*7) %>%
    group_by(occurence_week) %>%
    tally() %>%
    arrange(occurence_week) %>%
    ungroup() %>%
    mutate(delay_week = as.numeric((T-occurence_week)/7)) %>%
    left_join(delay_cdf, by = "delay_week") %>%
    mutate(log_rho = log(rho), t = 1:n())
  
  ##non-linear
  
  model_nb <- gam(
    n ~ s(t, bs = "tp") + offset(log_rho),
    family = nb(),
    method = "REML",
    data = regression_table
  )
  
  
  ###nowcast
  
  ### Just intensity
  pred_nowcast <- predict(
    model_nb,
    newdata = mutate(regression_table, log_rho = 0),
    type = "link",
    se.fit = TRUE,
    unconditional = TRUE
  )
  
  regression_table <- regression_table %>%
    mutate(
      fitted_reported = predict(
        model_nb,
        type = "response"
      ),
      nowcast = exp(pred_nowcast$fit),
      nowcast_rho = n / rho,
      nowcast_rho_low = n/rho_upper,
      nowcast_rho_up = n/rho_lower
    )
  
  #####data final
  
  data_final <- all_events %>%
    filter(occurence_week <= T_end, occurence_week >= min(regression_table$occurence_week)) %>%
    group_by(occurence_week) %>%
    summarise(true_n = n( )) %>%
    left_join(regression_table, by = "occurence_week")
  
  
  p <- ggplot(data_final, aes(x = occurence_week)) +
    geom_line(
      aes(y = true_n, color = "True"),
      linewidth = 1
    ) +
    geom_line(
      aes(y = nowcast, color = "Intensity"),
      linewidth = 1
    ) +
    geom_line(
      aes(y = n, color = "Observed"),
      linewidth = 1
    ) +
    geom_line(
      aes(y = nowcast_rho, color = "Nowcast Y/rho"),
      linewidth = 1
    )+
    geom_ribbon(
      aes(ymin = nowcast_rho_low, ymax = nowcast_rho_up),
      fill = "red",
      alpha = 0.2
    )+
    scale_color_manual(
      values = c(
        "Observed" = "black",
        "True" = "blue",
        "Intensity" = "green",
        "Nowcast Y/rho" = "red"
      )
    ) +
    labs(
      x = "Occurrence week",
      y = "Event count",
      color = NULL,
      title = paste0("Nowcasting ", country_name, " ", T)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom"
    )
  
  print(p)
  
}

dev.off()
