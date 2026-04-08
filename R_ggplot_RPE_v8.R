library(shiny)
library(dplyr)
library(readr)
library(ggplot2)
library(plotly)
library(viridis)
library(ggpointdensity)
library(forcats)

source("plotting_functions.R")

options(shiny.maxRequestSize = 8 * 1024^3)

mycolorv <- viridis(n = 11)
mycolorm <- magma(n = 11)
mycolorr <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

rename_if_present <- function(df, mapping) {
  keep <- mapping[names(mapping) %in% names(df)]
  if (length(keep) == 0) {
    return(df)
  }

  rename_args <- setNames(names(keep), unname(keep))
  dplyr::rename(df, !!!rename_args)
}

rename_channel_columns <- function(df, channels, aliases) {
  base_map <- list(
    DNA = c(
      int = "Intensity_IntegratedIntensity_DNA",
      intcor = "Intensity_IntegratedIntensity_dnacor",
      max = "Intensity_MaxIntensity_DNA",
      mean = "Intensity_MeanIntensity_DNA",
      meancor = "Intensity_MeanIntensity_dnacor",
      file = "FileName_DNA",
      path = "PathName_DNA"
    ),
    FITC = c(
      int = "Intensity_IntegratedIntensity_FITC",
      intcor = "Intensity_IntegratedIntensity_fitccor",
      mean = "Intensity_MeanIntensity_FITC",
      meancor = "Intensity_MeanIntensity_fitccor",
      file = "FileName_FITC",
      path = "PathName_FITC"
    ),
    mCherry = c(
      int = "Intensity_IntegratedIntensity_mCherry",
      intcor = "Intensity_IntegratedIntensity_mcherrycor",
      mean = "Intensity_MeanIntensity_mCherry",
      meancor = "Intensity_MeanIntensity_mcherrycor",
      file = "FileName_mCherry",
      path = "PathName_mCherry"
    ),
    Cy5 = c(
      int = "Intensity_IntegratedIntensity_Cy5",
      intcor = "Intensity_IntegratedIntensity_cy5cor",
      mean = "Intensity_MeanIntensity_Cy5",
      meancor = "Intensity_MeanIntensity_cy5cor",
      file = "FileName_Cy5",
      path = "PathName_Cy5"
    ),
    Cy3 = c(
      int = "Intensity_IntegratedIntensity_Cy3",
      intcor = "Intensity_IntegratedIntensity_cy3cor",
      mean = "Intensity_MeanIntensity_Cy3",
      meancor = "Intensity_MeanIntensity_cy3cor",
      file = "FileName_Cy3",
      path = "PathName_Cy3"
    )
  )

  for (ch in channels) {
    if (!ch %in% names(base_map)) {
      next
    }
    alias <- aliases[[ch]]
    for (nm in names(base_map[[ch]])) {
      src <- base_map[[ch]][[nm]]
      if (!src %in% names(df)) {
        next
      }
      suffix <- switch(
        nm,
        int = paste0("int_", alias),
        intcor = paste0("int_", alias, "cor"),
        max = paste0("max_", alias),
        mean = paste0("mean_", alias),
        meancor = paste0("mean_", alias, "cor"),
        file = paste0("file_", alias),
        path = paste0("path_", alias)
      )
      df <- dplyr::rename(df, !!suffix := all_of(src))
    }
  }

  df
}

process_pipeline <- function(nuc_data, plate_data, channels, aliases, var_names, var_order, cyto_data = NULL, cyto_var = NULL) {
  required_map <- c(
    Well = "Metadata_Well...6",
    col = "Metadata_col",
    row = "Metadata_row",
    scene = "Metadata_scene...11",
    area = "AreaShape_Area",
    convex = "AreaShape_ConvexArea",
    eccen = "AreaShape_Eccentricity",
    form = "AreaShape_FormFactor",
    length = "AreaShape_MajorAxisLength",
    min_feret = "AreaShape_MinFeretDiameter",
    min_length = "AreaShape_MinorAxisLength",
    image = "ImageNumber",
    nneigh = "Neighbors_NumberOfNeighbors_80",
    obj_id = "Number_Object_Number",
    x = "Location_Center_X",
    y = "Location_Center_Y",
    telophase = "Math_telophase"
  )

  nuc <- rename_if_present(nuc_data, required_map)
  nuc <- rename_channel_columns(nuc, channels, aliases)
  nuc <- mutate(nuc, i_id = paste(image, obj_id, sep = "_"))

  merged <- nuc
  if (!is.null(cyto_data) && !is.null(cyto_var) && nzchar(cyto_var)) {
    cyto_map <- c(
      Well = "Metadata_Well...6",
      col = "Metadata_col",
      row = "Metadata_row",
      scene = "Metadata_scene...11",
      image = "ImageNumber",
      obj_id = "Parent_Nuclei"
    )
    cyto <- rename_if_present(cyto_data, cyto_map)
    cyto <- rename_channel_columns(cyto, channels, aliases)
    cyto <- mutate(cyto, i_id = paste(image, obj_id, sep = "_"))

    cyto_col <- paste0("mean_", cyto_var, "cor")
    if (cyto_col %in% names(cyto) && cyto_col %in% names(nuc)) {
      cyto <- select(cyto, i_id, all_of(cyto_col))
      merged <- inner_join(nuc, cyto, by = "i_id", suffix = c("", "_cyto")) |>
        mutate(cn = .data[[paste0(cyto_col, "_cyto")]] / .data[[cyto_col]])
    }
  }

  dat <- inner_join(merged, plate_data, by = "Well") |>
    filter(form > 0.77)

  if (length(var_names) > 0 && length(var_order) >= length(var_names)) {
    dat <- FactorOrder(dat, var_names, var_order)
  }

  dat <- mutate(dat, scenetreat = interaction(Well, scene))

  if (length(var_names) > 1) {
    dat <- combinetreat(dat, var_names) |>
      mutate(treat = as.factor(treat))
  } else if (length(var_names) == 1) {
    dat <- mutate(dat, treat = .data[[var_names[1]]])
  } else {
    dat <- mutate(dat, treat = as.factor(Well))
  }

  if ("int_dnacor" %in% names(dat)) dat <- corany(dat, "int_dnacor", "scenetreat", "quart")
  if ("int_dna" %in% names(dat)) dat <- corany(dat, "int_dna", "scenetreat", "quart")

  for (ch in channels) {
    mean_cor <- paste0("mean_", aliases[[ch]], "cor")
    mean_raw <- paste0("mean_", aliases[[ch]])
    if (mean_cor %in% names(dat)) dat <- corany(dat, mean_cor, "scenetreat", "quant")
    if (mean_raw %in% names(dat)) dat <- corany(dat, mean_raw, "scenetreat", "quant")
  }

  if ("cor_int_dnacor" %in% names(dat)) {
    dat <- align(dat, "scenetreat", "int_dnacor", 0.2)
  }

  grouping_vars <- if (length(var_names) > 0 && all(var_names %in% names(dat))) {
    var_names
  } else {
    "treat"
  }

  hci.edu.count <- dat |>
    group_by(across(all_of(grouping_vars))) |>
    summarise(count = n(), .groups = "drop")

  min.count <- min(hci.edu.count$count)
  hci.edu.ds <- dat |>
    group_by(across(all_of(grouping_vars))) |>
    sample_n(min.count) |>
    ungroup()

  hci.edu.ds2 <- if (min.count > 3999) {
    dat |>
      group_by(across(all_of(grouping_vars))) |>
      sample_n(4000) |>
      ungroup()
  } else {
    hci.edu.ds
  }

  list(
    data = dat,
    count = hci.edu.count,
    min_count = min.count,
    downsampled = hci.edu.ds,
    downsampled_4k = hci.edu.ds2
  )
}

ui <- fluidPage(
  titlePanel("ImgCytometR Workflow App"),
  tabsetPanel(
    tabPanel(
      "1) Data Input",
      fluidRow(
        column(4,
          fileInput("nuc_file", "Nuclear CellProfiler CSV", accept = c(".csv")),
          fileInput("plate_file", "Plate Layout CSV", accept = c(".csv")),
          fileInput("cyto_file", "Cytoplasm CSV (optional)", accept = c(".csv")),
          checkboxInput("use_cyto", "Use cytoplasm join", value = FALSE),
          selectInput("channels", "Channels present", choices = c("DNA", "FITC", "mCherry", "Cy5", "Cy3"), multiple = TRUE,
                      selected = c("DNA", "FITC", "mCherry", "Cy5")),
          actionButton("run_pipeline", "Run preprocessing", class = "btn-primary")
        ),
        column(8,
          h4("Processed data preview"),
          tableOutput("processed_preview")
        )
      )
    ),
    tabPanel(
      "2) Naming",
      fluidRow(
        column(3,
          textInput("alias_DNA", "DNA alias", value = "dna"),
          textInput("alias_FITC", "FITC alias", value = "edu"),
          textInput("alias_mCherry", "mCherry alias", value = "cdt1"),
          textInput("alias_Cy5", "Cy5 alias", value = "cdt1_ha"),
          textInput("alias_Cy3", "Cy3 alias", value = "cy3"),
          textInput("cyto_var", "Cytoplasm variable alias", value = "prb1")
        ),
        column(4,
          textInput("var_names", "Treatment columns (comma-separated)", value = "treat_1,treat_2,treat_3"),
          textInput("var_order", "Order columns (comma-separated)", value = "order_1,order_2,order_3"),
          helpText("These drive combined treatment naming and ordering.")
        ),
        column(5,
          h4("Current settings"),
          verbatimTextOutput("naming_settings")
        )
      )
    ),
    tabPanel(
      "3) Gating",
      fluidRow(
        column(3,
          numericInput("x1", "x1", value = 25),
          numericInput("x2", "x2", value = 98),
          numericInput("x5", "x5", value = 81),
          numericInput("x6", "x6", value = 120),
          numericInput("y2", "y2", value = 0.026)
        ),
        column(9,
          plotlyOutput("gate_plot", height = "600px"),
          tableOutput("phase_summary")
        )
      )
    ),
    tabPanel(
      "4) Existing Shiny Apps",
      p("These buttons call the previously defined apps in plotting_functions.R."),
      actionButton("launch_2d", "Launch shiny2D()"),
      actionButton("launch_violin", "Launch shinyViolin()"),
      actionButton("launch_percent", "Launch shinyPercent()")
    )
  )
)

server <- function(input, output, session) {
  aliases <- reactive({
    list(
      DNA = input$alias_DNA,
      FITC = input$alias_FITC,
      mCherry = input$alias_mCherry,
      Cy5 = input$alias_Cy5,
      Cy3 = input$alias_Cy3
    )
  })

  parsed_var_names <- reactive({
    x <- trimws(unlist(strsplit(input$var_names, ",")))
    x[x != ""]
  })

  parsed_var_order <- reactive({
    x <- trimws(unlist(strsplit(input$var_order, ",")))
    x[x != ""]
  })

  processed_bundle <- eventReactive(input$run_pipeline, {
    req(input$nuc_file, input$plate_file)

    nuc <- readr::read_csv(input$nuc_file$datapath, show_col_types = FALSE)
    plate <- readr::read_csv(input$plate_file$datapath, show_col_types = FALSE)
    cyto <- NULL
    if (isTRUE(input$use_cyto) && !is.null(input$cyto_file)) {
      cyto <- readr::read_csv(input$cyto_file$datapath, show_col_types = FALSE)
    }

    process_pipeline(
      nuc_data = nuc,
      plate_data = plate,
      channels = input$channels,
      aliases = aliases(),
      var_names = parsed_var_names(),
      var_order = parsed_var_order(),
      cyto_data = cyto,
      cyto_var = input$cyto_var
    )
  })

  gated <- reactive({
    dat <- processed_bundle()$downsampled
    req("cor_int_dnacor" %in% names(dat), "cor_mean_educor" %in% names(dat))

    x1 <- input$x1
    x2 <- input$x2
    x3 <- x2 + (x2 - x1)
    x5 <- input$x5
    x6 <- input$x6
    y1 <- min(dat$cor_mean_educor, na.rm = TRUE)
    y2 <- input$y2
    y4 <- max(dat$cor_mean_educor, na.rm = TRUE)

    gated.edu.f <- gating(dat, "cor_mean_educor", "cor_int_dnacor")
    gated.edu.f$phase2 <- factor(gated.edu.f$phase2, levels = c("G1", "S", "G2"))
    gated.edu.f
  })

  output$processed_preview <- renderTable({
    head(processed_bundle()$data, 10)
  })

  output$naming_settings <- renderPrint({
    list(
      channels = input$channels,
      aliases = aliases(),
      treatment_columns = parsed_var_names(),
      order_columns = parsed_var_order(),
      cyto_var = input$cyto_var
    )
  })

  output$gate_plot <- renderPlotly({
    dat <- processed_bundle()$downsampled
    req("cor_int_dnacor" %in% names(dat), "cor_mean_educor" %in% names(dat))

    x3 <- input$x2 + (input$x2 - input$x1)
    y1 <- min(dat$cor_mean_educor, na.rm = TRUE)
    y4 <- max(dat$cor_mean_educor, na.rm = TRUE)

    p <- ggplot(dat, aes(x = cor_int_dnacor, y = cor_mean_educor)) +
      geom_pointdensity(size = 0.3, adjust = 0.05) +
      scale_y_continuous(trans = "log10") +
      scale_color_gradientn(colours = mycolorv) +
      annotate("rect", xmin = input$x1, xmax = input$x2, ymin = y1, ymax = input$y2, alpha = 0, color = "red", linetype = "dotted") +
      annotate("rect", xmin = input$x1, xmax = input$x5, ymin = input$y2, ymax = y4, alpha = 0, color = "red", linetype = "dotted") +
      annotate("rect", xmin = input$x5, xmax = input$x6, ymin = input$y2, ymax = y4, alpha = 0, color = "red", linetype = "dotted") +
      annotate("rect", xmin = input$x6, xmax = x3, ymin = input$y2, ymax = y4, alpha = 0, color = "red", linetype = "dotted") +
      annotate("rect", xmin = input$x2, xmax = x3, ymin = y1, ymax = input$y2, alpha = 0, color = "red", linetype = "dotted") +
      theme_classic()

    ggplotly(p)
  })

  output$phase_summary <- renderTable({
    g <- gated()
    g |>
      group_by(treat, phase2) |>
      summarise(count = n(), .groups = "drop_last") |>
      mutate(per = prop.table(count) * 100) |>
      ungroup()
  })

  observeEvent(input$launch_2d, {
    shiny::runApp(shiny2D())
  })

  observeEvent(input$launch_violin, {
    shiny::runApp(shinyViolin())
  })

  observeEvent(input$launch_percent, {
    shiny::runApp(shinyPercent())
  })
}

shinyApp(ui, server)
