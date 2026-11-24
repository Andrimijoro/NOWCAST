library(mgcv)
library(dplyr)
library(ggplot2)

set.seed(123)

B <- 2000L

risk_data <- risk_data %>%
  ungroup() %>%
  arrange(delay_week)

# ============================================================
# 1. ORIGINAL DELAY MODEL
# ============================================================

delay_model <- mgcv::gam(
  cbind(hits, risk - hits) ~
    s(delay_week, bs = "tp"),
  data = risk_data,
  family = binomial(link = "cloglog"),
  method = "REML"
)

delay_grid <- risk_data$delay_week
K <- length(delay_grid)

# Initial population at risk
N0 <- risk_data$risk[1]

# Check that risk sets are consistent
stopifnot(
  all(diff(delay_grid) == 1),
  min(delay_grid) == 0,
  all(
    risk_data$risk[-1] ==
      head(risk_data$risk - risk_data$hits, -1)
  )
)

# Estimated hazard
hazard_hat <- predict(
  delay_model,
  newdata = risk_data,
  type = "response"
)

# Reporting completeness
rho_hat <- 1 - cumprod(1 - hazard_hat)

# ============================================================
# 2. PREPARE NOWCAST DATA
# ============================================================

nowcast_table <- regression_table %>%
  mutate(
    age = floor(
      as.numeric(T - as.Date(occurence_week)) / 7
    )
  ) %>%
  filter(age >= 0) %>%
  mutate(
    delay_index = match(
      pmin(age, max(delay_grid)),
      delay_grid
    ),
    rho = rho_hat[delay_index],
    nowcast = n / rho
  )

# ============================================================
# 3. BOOTSTRAP
# ============================================================

rho_boot <- matrix(
  NA_real_,
  nrow = nrow(nowcast_table),
  ncol = B
)

nowcast_boot <- rho_boot

for (b in seq_len(B)) {
  
  hits_b <- numeric(K)
  risk_b <- numeric(K)
  
  risk_b[1] <- N0
  
  # Simulate reporting process sequentially
  for (j in seq_len(K)) {
    
    hits_b[j] <- rbinom(
      n = 1,
      size = risk_b[j],
      prob = hazard_hat[j]
    )
    
    if (j < K) {
      risk_b[j + 1] <- risk_b[j] - hits_b[j]
    }
  }
  
  risk_data_b <- risk_data %>%
    mutate(
      hits = hits_b,
      risk = risk_b
    ) %>%
    filter(risk > 0)
  
  # Refit delay model
  model_b <- tryCatch(
    mgcv::gam(
      cbind(hits, risk - hits) ~
        s(delay_week, bs = "tp"),
      data = risk_data_b,
      family = binomial(link = "cloglog"),
      method = "REML"
    ),
    error = function(e) NULL
  )
  
  if (is.null(model_b)) next
  
  # Bootstrap hazard
  hazard_b <- predict(
    model_b,
    newdata = data.frame(
      delay_week = delay_grid
    ),
    type = "response"
  )
  
  # Bootstrap completeness
  completeness_b <- 1 - cumprod(1 - hazard_b)
  
  rho_b <- completeness_b[
    nowcast_table$delay_index
  ]
  
  if (any(!is.finite(rho_b) | rho_b <= 0)) next
  
  # Save bootstrap estimates
  rho_boot[, b] <- rho_b
  
  nowcast_boot[, b] <- nowcast_table$n / rho_b
}

# Retain successful bootstrap iterations
valid <- colSums(is.finite(nowcast_boot)) ==
  nrow(nowcast_table)

rho_boot <- rho_boot[, valid, drop = FALSE]
nowcast_boot <- nowcast_boot[, valid, drop = FALSE]

stopifnot(ncol(nowcast_boot) > 1)

# ============================================================
# 4. BOOTSTRAP CONFIDENCE INTERVALS
# ============================================================

nowcast_table <- nowcast_table %>%
  mutate(
    lower = apply(
      nowcast_boot, 1, quantile, probs = 0.025
    ),
    upper = apply(
      nowcast_boot, 1, quantile, probs = 0.975
    ),
    bootstrap_se = apply(
      nowcast_boot, 1, sd
    ),
    bootstrap_mean = rowMeans(nowcast_boot)
  )

ggplot(
  nowcast_table,
  aes(x = occurence_week)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "red",
    alpha = 0.15
  ) +
  geom_line(
    aes(y = nowcast),
    color = "red",
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = n),
    color = "black"
  ) +
  labs(
    x = "Occurrence week",
    y = "Event count",
    title = "Nowcast with 95% bootstrap confidence intervals"
  ) +
  theme_minimal()


nowcast_table <- nowcast_table %>%
  mutate(
    rho_se = apply(rho_boot, 1, sd),
    
    var_delay = (n^2 / rho^4) * rho_se^2,
    
    var_reporting = nowcast * (1 - rho) / rho,
    
    var_total = var_delay + var_reporting,
    
    se_total = sqrt(var_total),
    
    lower_total = pmax(n, nowcast - 1.96 * se_total),
    
    upper_total = nowcast + 1.96 * se_total
  )


nowcast_table <- nowcast_table %>%
  mutate(
    nowcast = n / rho,
    
    lower = ifelse(
      n == 0, 0,
      n + qnbinom(
        0.025,
        size = pmax(n, 1),
        prob = rho
      )
    ),
    
    upper = n + qnbinom(
      0.975,
      size = n + 1,
      prob = rho
    )
  )
