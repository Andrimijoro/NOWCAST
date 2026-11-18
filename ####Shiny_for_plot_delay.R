library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# Load data
data <- readRDS("base_df.rds")

data <- data %>%
  mutate(
    timestamp_date = as.Date(
      as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC")
    ),
    delay = floor(as.numeric(timestamp_date - event_date) / 7),
    occurence_week = floor_date(
      event_date, unit = "week", week_start = 1
    ),
    report_week = floor_date(
      timestamp_date, unit = "week", week_start = 1
    )
  ) %>%
  mutate(delay = ifelse(delay > 20, 20, delay)) %>%
  ungroup()

ui <- fluidPage(
  
  titlePanel("Reporting Delays in Africa"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "country",
        label = "Select country",
        choices = c(
          "All countries",
          sort(unique(data$country))
        ),
        selected = "All countries"
      ),
      
      downloadButton(
        outputId = "download_plot",
        label = "Download plot",
        icon = icon("download")
      ),
      
      width = 3
    ),
    
    mainPanel(
      plotOutput("delay_plot", height = "550px"),
      width = 9
    )
  )
)

server <- function(input, output, session) {
  
  filtered_data <- reactive({
    
    df <- data
    
    if (input$country != "All countries") {
      df <- df %>%
        filter(country == input$country)
    }
    
    df
  })
  
  delay_summary <- reactive({
    
    filtered_data() %>%
      filter(
        occurence_week <=
          max(occurence_week, na.rm = TRUE) - 20 * 7
      ) %>%
      group_by(occurence_week) %>%
      summarise(
        n = n(),
        median = median(delay, na.rm = TRUE),
        q25 = quantile(delay, 0.25, na.rm = TRUE),
        q75 = quantile(delay, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  # Reusable plot
  delay_plot <- reactive({
    
    df_counts <- delay_summary()
    
    # Scaling factor for secondary y-axis
    scale_factor <- max(df_counts$n, na.rm = TRUE) /
      max(1, df_counts$q75, na.rm = TRUE)
    
    df1 <- df_counts %>%
      pivot_longer(
        cols = c(median, q25, q75),
        names_to = "Statistic",
        values_to = "Delay"
      ) %>%
      mutate(
        Statistic = factor(
          Statistic,
          levels = c("median", "q75", "q25"),
          labels = c(
            "Median",
            "75th percentile",
            "25th percentile"
          )
        )
      )
    
    ggplot() +
      
      # Unsmoothed reporting-delay quantiles
      geom_line(
        data = df1,
        aes(
          x = occurence_week,
          y = Delay,
          color = Statistic,
          linetype = Statistic
        ),
        linewidth = 0.9
      ) +
      
      # Smoothed weekly event counts
      geom_smooth(
        data = df_counts,
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
        ),
        breaks = c(
          "Median",
          "75th percentile",
          "25th percentile",
          "Number of events"
        )
      ) +
      
      scale_linetype_manual(
        values = c(
          "Median" = "solid",
          "75th percentile" = "dashed",
          "25th percentile" = "dashed",
          "Number of events" = "solid"
        ),
        breaks = c(
          "Median",
          "75th percentile",
          "25th percentile",
          "Number of events"
        )
      ) +
      
      # Dual y-axis
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
        axis.title.y.right = element_text(
          color = "#555555"
        ),
        legend.key.width = grid::unit(1.5, "cm")
      )
  })
  
  # Display plot
  output$delay_plot <- renderPlot({
    delay_plot()
  })
  
  # Download PDF
  output$download_plot <- downloadHandler(
    
    filename = function() {
      
      country_name <- if (input$country == "All countries") {
        "Africa"
      } else {
        input$country
      }
      
      country_name <- gsub(
        "[^A-Za-z0-9_-]", "_", country_name
      )
      
      paste0("reporting_delay_", country_name, ".pdf")
    },
    
    content = function(file) {
      
      ggsave(
        filename = file,
        plot = delay_plot(),
        device = "pdf",
        width = 10,
        height = 6,
        units = "in",
        bg = "white"
      )
    },
    
    contentType = "application/pdf"
  )
}

shinyApp(ui = ui, server = server)
