FactorOrder <- function(df, varNames, varOrder){
  for(i in 1:length(varNames)) {
    df <- df %>% arrange(.data[[varOrder[i]]]) %>% 
      dplyr::mutate("{varNames[i]}" := factor(.data[[varNames[i]]],unique(.data[[varNames[i]]])))
  }
  return(df)
}

combinetreat <- function(df, varNames){
  df <- df %>% mutate(treat = paste(.data[[varNames[1]]],.data[[varNames[2]]],sep="+"))
  if(length(varNames) > 2){
    for(i in 3:length(varNames)){
      df <- df %>% mutate(treat = paste(treat,.data[[varNames[i]]]))
    }
  }
  return(df)
}

corany <- function(a,b,c,d){
  if(d == "mean"){
    corfac <- a  %>% group_by(.data[[c]]) %>% summarise(mean = mean(.data[[b]])) %>%
      mutate(corfac = mean[1]/mean) %>% dplyr::select(-mean) 
    
    cor <- left_join(x = a, y = corfac) %>% mutate("cor_{b}" := (.data[[b]]*corfac)) %>% dplyr::select(-corfac)                 
    return(cor)
  }
  if(d == "min"){
    corfac <- a  %>% group_by(.data[[c]]) %>% summarise(mean = min(.data[[b]])) %>%
      mutate(corfac = mean[1]/mean) %>% dplyr::select(-mean)  
    
    cor <- left_join(x = a, y = corfac) %>% mutate("cor_{b}" := (.data[[b]]*corfac)) %>% dplyr::select(-corfac)                 
    return(cor)
  }
  if(d == "quant"){
    corfac <- a  %>% group_by(.data[[c]]) %>% summarise(mean = quantile(.data[[b]], 0.01)) %>%
      mutate(corfac = mean[1]/mean) %>% dplyr::select(-mean)
    
    cor <- left_join(x = a, y = corfac) %>% mutate("cor_{b}" := (.data[[b]]*corfac)) %>% dplyr::select(-corfac)                 
    return(cor)
  }
  if(d == "quart"){
    corfac <- a  %>% group_by(.data[[c]]) %>% summarise(mean = quantile(.data[[b]], 0.1)) %>%
      mutate(corfac = mean[1]/mean) %>% dplyr::select(-mean)
    
    cor <- left_join(x = a, y = corfac) %>% mutate("cor_{b}" := (.data[[b]]*corfac)) %>% dplyr::select(-corfac)                 
    return(cor)
  }
  else{return(NA)}
}

align <- function(a,b,c,ad){          
  df <- droplevels(a)
  df[[b]] <- as.factor(df[[b]])
  e <- list()
  for(i in 1:length(levels(df[[b]]))){
    f <- df %>% filter(.data[[b]] == levels(df[[b]])[i]) 
    if (nrow(f) == 1){
      e[[b]][i] <- levels(df[[b]])[i]
      e$g1peak[i] <- f[[c]]
      next}
    d <- density(f[[c]], adjust = ad)
    max <- d$x[d$y==max(d$y)]
    
    
    e[[b]][i] <- levels(df[[b]])[i]
    e$g1peak[i] <- max[1]
  }
  
  e <- as.data.frame(e)
  corfac <- e %>% mutate(corfac = g1peak[1]/g1peak) %>% dplyr::select(-g1peak) 
  cor <- left_join(x = a, y = corfac) %>% mutate("a{c}" := (.data[[c]]*corfac)) %>% dplyr::select(-corfac) 
  
  return(cor)
}

gating <- function(df, y_var,x_var){
  
  hci.edu.tib <- df %>% as_tibble
  
  gated.hci.edu <- hci.edu.tib %>% mutate(phase = ifelse(between(.data[[y_var]],y1,y2) & between(.data[[x_var]] ,x1,x2), "G1",
                                                         ifelse(between(.data[[y_var]],y2,y4) & between(.data[[x_var]] ,x1,x2), "ES",
                                                                ifelse(between(.data[[y_var]],y2,y4) & between(.data[[x_var]] ,x2,x3), "LS",
                                                                       ifelse(between(.data[[y_var]],y1,y2) & between(.data[[x_var]] ,x2,x3), "G2", "OT"))))) 
  
  tgated.hci.edu <- gated.hci.edu %>% mutate(phase2 = ifelse(between(.data[[y_var]],y1,y2) & between(.data[[x_var]] ,x1,x2), "G1",
                                                             ifelse(between(.data[[y_var]],y2,y4) & between(.data[[x_var]] ,x1,x3), "S",
                                                                    ifelse(between(.data[[y_var]],y1,y2) & between(.data[[x_var]] ,x2,x3), "G2", "OT"))))
  
  fivegated.hci.edu <- tgated.hci.edu %>% mutate(phase5 = ifelse(between(.data[[y_var]],y1,y2) & between(.data[[x_var]] ,x1,x2), "G1",
                                                                 ifelse(between(.data[[y_var]],y2,y4) & between(.data[[x_var]] ,x1,x5), "ES",
                                                                        ifelse(between(.data[[y_var]],y2,y4) & between(.data[[x_var]] ,x5,x6), "MS",
                                                                               ifelse(between(.data[[y_var]],y2,y4) & between(.data[[x_var]] ,x6,x3), "LS",
                                                                                      ifelse(between(.data[[y_var]],y1,y2) & between(.data[[x_var]] ,x2,x3), "G2", "OT")))))) 
  
  fivegated.hci.edu %>%filter (phase != "OT") %>%filter(phase2 != "OT") %>% filter(phase5 != "OT")
  
}

telo_id <- function(telo){
telo_cord <- telo %>% dplyr::select(c("x","y", "i_id", "image"))
matches <- list()
x <- 1
for (a in 1:nrow(telo)){
  filt_telo_cord <- telo_cord %>% filter(image == telo$image[a]) %>% filter(i_id != telo$i_id[a]) 
  if (nrow(filt_telo_cord) == 0){next}
  else{    
    for (b in 1:nrow(filt_telo_cord)){
      if (between(telo$x[a], (filt_telo_cord$x[b]-35), (filt_telo_cord$x[b]+35)) & between(telo$y[a], (filt_telo_cord$y[b]-35), (filt_telo_cord$y[b]+35))){
        matches$parent[x] <- telo$i_id[a]
        matches$x[x] <- telo$x[a]
        matches$y[x] <- telo$y[a]
        matches$pair[x] <- filt_telo_cord$i_id[b]
        matches$x_pair[x] <- filt_telo_cord$x[b]
        matches$y_pair[x] <- filt_telo_cord$y[b]
        matches$match[x] <- T
      } 
      else {
        matches$parent[x] <- telo$i_id[a]
        matches$x[x] <- telo$x[a]
        matches$y[x] <- telo$y[a]
        matches$pair[x] <- filt_telo_cord$i_id[b]
        matches$x_pair[x] <- filt_telo_cord$x[b]
        matches$y_pair[x] <- filt_telo_cord$y[b]
        matches$match[x] <- F
      }
      x <- x+1
    }
  }
  remove(filt_telo_cord)
}
return(matches)
}

foldtest <- function(a,b,c,d,e){
  a %>% mutate(
    "fold_{d}" := .data[[b]]/c,
    "{d}_test" := ifelse(.data[[paste0("fold_",d)]] > e, "Pos+", "Neg-"))
}

percents <- function(df, groups){ 
  loaded <- df %>% filter() %>%
    group_by_at(groups) %>% summarise(count = n()) %>% mutate(per= prop.table(count) * 100)
  return(loaded)
}


shiny2D <- function(){
shinyApp(
  ui <- fluidPage(
    tabsetPanel(tabPanel("Plots",
                         fluidRow(column(12,
                                         plotOutput("plot", inline = TRUE))),
                         fluidRow(column(2,
                                         actionButton("graph_click", "Make Graph"),
                                         selectInput("dataset", "Select Data Set", c( "Full" = "hci.edu.c",
                                                                                      "Large Downsample" = "hci.edu.ds",
                                                                                      "Small Downsample" = "hci.edu.ds2",
                                                                                      "Gated" = "gated.edu.f"
                                         )),
                                         varSelectInput("x_var", "X_Axis Variable", dplyr::select(hci.edu.c,where(is.numeric)),
                                                        selected = "cor_int_dnacor"),
                                         varSelectInput("y_var", "Y_Axis Variable", dplyr::select(hci.edu.c,where(is.numeric)),
                                                        selected = "cor_mean_educor"),
                                         radioButtons("faceting", "Facet Type", c("None", "Wrap", "Grid")),
                                         selectInput("colFacet", "Column/Wrap Facet Var:", c("None", "treat", varNames, "phase", "phase2", "phase5")),
                                         selectInput("colRow", "Row Facet Var:", c("None", "treat", varNames, "phase", "phase2", "phase5"))
                         ),
                         
                         column(3,
                                sliderInput("height", "height", min = 100, max = 2000, value = 500),
                                sliderInput("width", "width", min = 100, max = 2000, value = 1500),
                                numericInput("res", "Scaling (click Make Graph after change)", value = 64),
                                numericInput("dot_size", "Point Size", value = 0.3),
                                numericInput("gradient_size", "Change density gradient", value = 0.01),
                                numericInput("text_globe", "Global text size", value = 28),
                                numericInput("text_axis", "Axis text size", value = 24),
                                numericInput("text_leg", "Legend text size", value = 20),
                         ),
                         column(2,
                                radioButtons("filtering", "Filter groups?", c("Yes", "No"), selected = "No"),
                                checkboxGroupInput("levels", "Groups to include:", c(levels(as.factor(hci.edu.c$treat)))),
                                radioButtons("phase_filtering", "Filter phases?", c("Yes", "No"), selected = "No"),
                                checkboxGroupInput("phases", "Phases to include:", if(exists("gated.edu.f")){c(levels(as.factor(gated.edu.f$phase5)))} else{c("No phases")})
                         ),
                         column(3,
                                selectInput("tran_x", "X-Axis Transformation", c("Linear" = "identity",
                                                                                 "Log10" = "log10",
                                                                                 "pseudolog" = "pseudo_log"),
                                            selected = "identity"),
                                numericInput("sigma_x", "X-axis Sigma (For pseudolog only)", value = NULL),
                                selectInput("tran_y", "Y-Axis Transformation", c("Linear" = "identity",
                                                                                 "Log10" = "log10",
                                                                                 "pseudolog" = "pseudo_log"),
                                            selected = "log10"),
                                numericInput("sigma_y", "Y-axis Sigma (For pseudolog only)", value = NULL),
                                textInput("title", "Title"),
                                textInput("x_lab", "X-label"),
                                textInput("y_lab", "Y-label"),
                                selectInput("twocolor", "Select Color Map", c("Magma" = "mycolorm",
                                                                              "Rainbow" = "mycolorr",
                                                                              "Viridis" = "mycolorv"),
                                            selected = "mycolorv")),
                         column(2,
                                numericInput("lower_x", "Lower X Limit", value = NULL),
                                numericInput("upper_x", "Upper X Limit", value = NULL),
                                numericInput("lower_y", "Lower Y Limit", value = NULL),
                                numericInput("upper_y", "Upper Y Limit", value = NULL),
                                numericInput("aspect", "Aspect Ratio (No Number for free axis)", value = NULL)
                         )),
                         fluidRow(column(3,
                                         textInput("graph_file", "Graph File Name"),
                                         verbatimTextOutput("code")),                               ####Debug
                                  column(2,
                                         downloadButton("graph", "Save Graph (Use .pdf extensions):")))
                         
    )
    )
  )
  ,
  
  server <- function(input, output, session) {
    gg <- reactive(  
      get(input$dataset) %>% 
        {if(input$filtering == "Yes") filter(.,treat %in% c(input$levels)) else .}  %>%
        {if(input$phase_filtering == "Yes") filter(.,phase5 %in% c(input$phases)) else .}  %>% 
        ggplot(aes(x=!!input$x_var, y=!!input$y_var))+
        geom_pointdensity(size = input$dot_size, 
                          adjust = input$gradient_size) +
        scale_y_continuous(transform = if(input$tran_y == "pseudo_log"){scales::pseudo_log_trans(sigma = input$sigma_y)}
                           else{input$tran_y},
                           limits = c(input$lower_y, 
                                      input$upper_y))+
        scale_x_continuous(transform = if(input$tran_x == "pseudo_log"){scales::pseudo_log_trans(sigma = input$sigma_x)}
                           else{input$tran_x},
                           limits = c(input$lower_x, 
                                      input$upper_x))+
        scale_color_gradientn(colours = get(input$twocolor))+
        labs(title = input$title,
             x= input$x_lab,
             y=input$y_lab, 
             colour = "Density")+
        theme_classic()+
        theme(text = element_text(size=input$text_globe),
              legend.title = element_text(size = input$text_leg),
              axis.text = element_text(size=input$text_axis),
              plot.title = element_text(hjust = 0.5),
              axis.line = element_blank(),
              strip.background = element_blank(),
              panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
        {if(is.na(input$aspect)) theme(aspect.ratio=NULL) else theme(aspect.ratio = input$aspect)}+
        {if(input$faceting == "Grid")
          facet_grid(cols = eval(if(input$colFacet == "None") NULL else vars(get(input$colFacet))),
                     rows = eval(if(input$colRow == "None") NULL else vars(get(input$colRow))))}+
        {if(input$faceting == "Wrap")
          facet_wrap(eval(if(input$colFacet == "None") NULL else vars(get(input$colFacet))))} 
    )
    
    observeEvent(input$graph_click,{
      output$plot <- renderPlot({
        gg()
      },
      width = function() input$width,
      height = function() input$height,
      res = input$res
      )
    })
    output$graph <- downloadHandler(filename = function(){
      input$graph_file}, 
      content =  function(file) {
        ggsave(file, gg(), device = "pdf",
               width = input$width,
               height = input$height,
               units = c("px"),
               dpi = input$res) 
        #pdf(file=file)
        #plot(gg()) 
        #dev.off()
      }
    )
  }
)
}

##Violin Shiny#
shinyViolin <- function(){ 
shinyApp(
    ui <- fluidPage(
      tabsetPanel(tabPanel("Plots",
                           fluidRow(column(12,
                                           plotOutput("plot", inline = TRUE))),
                           fluidRow(column(2,
                                           actionButton("graph_click", "Make Graph"),
                                           selectInput("dataset", "Select Data Set", c( "Full" = "hci.edu.c",
                                                                                        "Large Downsample" = "hci.edu.ds",
                                                                                        "Small Downsample" = "hci.edu.ds2",
                                                                                        "Gated" = "gated.edu.f",
                                                                                        "Telophases" = "v.unpair"
                                           )),
                            varSelectInput("x_var", "X_Axis Variable", if(exists("gated.edu.f")){dplyr::select(gated.edu.f, where(is.factor))} else{dplyr::select(hci.edu.c, where(is.factor))}, selected = "treat"),
                            varSelectInput("y_var", "Y_Axis Variable", if(exists("gated.edu.f")){gated.edu.f} else{hci.edu.c}, selected = "cor_mean_educor"),
                                           radioButtons("faceting", "Facet Type", c("None", "Wrap", "Grid")),
                                           selectInput("colFacet", "Column/Wrap Facet Var:", c("None", "treat", varNames, "phase", "phase2", "phase5")),
                                           selectInput("colRow", "Row Facet Var:", c("None", "treat", varNames, "phase", "phase2", "phase5"))
                           ),
                           
                           column(3,
                                  sliderInput("height", "height", min = 100, max = 2000, value = 600),
                                  sliderInput("width", "width", min = 100, max = 2000, value = 600),
                                  numericInput("res", "Scaling (click Make Graph after change)", value = 64),
                                  numericInput("vio_alpha", "Violin Transparency", value = 0.2),
                                  numericInput("dot_size", "Point Size", value = 0.1),
                                  numericInput("dot_alpha", "Dot Transparency", value = 0.1)

                           ),
                           column(2,
                                  radioButtons("filtering", "Filter groups?", c("Yes", "No"), selected = "No"),
                                  checkboxGroupInput("levels", "Groups to include:", c(levels(as.factor(hci.edu.c$treat)))),
                                  radioButtons("phase_filtering", "Filter phases?", c("Yes", "No"), selected = "No"),
                                  checkboxGroupInput("phases", "Phases to include:", if(exists("gated.edu.f")){c(levels(as.factor(gated.edu.f$phase5)))} else{c("No phases")})
                           ),
                           column(3,
                                  selectInput("tran_y", "Y-Axis Transformation", c("Linear" = "identity",
                                                                                   "Log10" = "log10",
                                                                                   "pseudolog" = "pseudo_log"),
                                              selected = "log10"),
                                  numericInput("sigma_y", "Y-axis Sigma (For pseudolog only)", value = NULL),
                                  textInput("title", "Title"),
                                  textInput("x_lab", "X-label"),
                                  textInput("y_lab", "Y-label"),
                                  textInput("fill", "Select Color", value = "deepskyblue4"),
                                  radioButtons("telo_filtering", "Late Telohase filter?", c("Yes", "No"), selected = "No")),
                           column(2,
                                  numericInput("lower_y", "Lower Y Limit", value = NULL),
                                  numericInput("upper_y", "Upper Y Limit", value = NULL),
                                  numericInput("aspect", "Aspect Ratio (No Number for free axis)", value = 2),
                                  numericInput("text_globe", "Global text size", value = 28),
                                  numericInput("text_axis", "Axis text size", value = 24),
                                  numericInput("angle", "X-text Angle", value = 45),
                                  numericInput("hjust", "X-text alignment", value = 1)
                           )),
                           fluidRow(column(4,
                                           textInput("graph_file", "Graph File Name"),
                                           verbatimTextOutput("code")),                               ####Debug
                                    column(3,
                                           downloadButton("graph", "Save Graph (Use .pdf extensions):"))),
                           fluidRow(column(12, 
                           p("Color names can be found here: https://sape.inf.usi.ch/sites/default/files/ggplot2-colour-names.png")))
                           
      )
      )
    )
    ,
    
    server <- function(input, output, session) {
      gg <- reactive(  
        get(input$dataset) %>% 
          {if(input$filtering == "Yes") filter(.,treat %in% c(input$levels)) else .}  %>%
          {if(input$phase_filtering == "Yes") filter(.,phase5 %in% c(input$phases)) else .}  %>% 
          {if(input$telo_filtering == "Yes") filter(.,telophase < late_telo) else .}  %>% 
          ggplot(aes(x=!!input$x_var, y=!!input$y_var))+
          geom_violin(linewidth = 1, scale = "width", fill = input$fill, alpha = input$vio_alpha, draw_quantiles = c(0.5))+
          geom_quasirandom(size = input$dot_size, alpha= input$dot_alpha)+
          scale_y_continuous(transform = if(input$tran_y == "pseudo_log"){scales::pseudo_log_trans(sigma = input$sigma_y)}
                             else{input$tran_y},
                             limits = c(input$lower_y, 
                                        input$upper_y))+
          labs(title = if(input$telo_filtering == "Yes"){"Late Telophase"} else{input$title},
               x= input$x_lab,
               y=input$y_lab)+
          theme_classic()+
          theme(text = element_text(size=input$text_globe),
                legend.title = element_text(size = input$text_leg),
                axis.text.y = element_text(size=input$text_axis),
                axis.text.x = element_text(size=input$text_axis, angle = input$angle, hjust=input$hjust),
                plot.title = element_text(hjust = 0.5),
                axis.line = element_blank(),
                strip.background = element_blank(),
                panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
          {if(is.na(input$aspect)) theme(aspect.ratio=NULL) else theme(aspect.ratio = input$aspect)}+
          {if(input$faceting == "Grid")
            facet_grid(cols = eval(if(input$colFacet == "None") NULL else vars(get(input$colFacet))),
                       rows = eval(if(input$colRow == "None") NULL else vars(get(input$colRow))))}+
          {if(input$faceting == "Wrap")
            facet_wrap(eval(if(input$colFacet == "None") NULL else vars(get(input$colFacet))))} 
      )
      
      observeEvent(input$graph_click,{
        output$plot <- renderPlot({
          gg()
        },
        width = function() input$width,
        height = function() input$height,
        res = input$res
        )
      })
      output$graph <- downloadHandler(filename = function(){
        input$graph_file}, 
        content =  function(file) {
          ggsave(file, gg(), device = "pdf",
                 width = input$width,
                 height = input$height,
                 units = c("px"),
                 dpi = input$res) 
          #pdf(file=file)
          #plot(gg()) 
          #dev.off()
        }
      )
    }
  )
}

#Percent Shiny#
shinyPercent <- function(){
shinyApp(
    ui <- fluidPage(
      tabsetPanel(tabPanel("Set Gate",
                           fluidRow(column(12,
                                           plotOutput("plot", inline = TRUE))),
                           fluidRow(column(2,
                                           actionButton("graph_click", "Make Graph"),
                                           selectInput("dataset", "Select Data Set", c( "Gated" = "gated.edu.f",
                                                                                        "Telophases" = "v.unpair",
                                                                                        "Late Telophase" = "late_telophase"
                                           )),
                                           varSelectInput("x_var", "X_Axis Variable", dplyr::select(hci.edu.c,where(is.numeric)),
                                                          selected = "cor_int_dnacor"),
                                           varSelectInput("y_var", "Y_Axis Variable", dplyr::select(hci.edu.c,where(is.numeric)),
                                                          selected = "cor_mean_educor"),
                                           radioButtons("faceting", "Facet Type", c("None", "Wrap", "Grid")),
                                           selectInput("colFacet", "Column/Wrap Facet Var:", c("None", "treat", varNames, "phase", "phase2", "phase5")),
                                           selectInput("colRow", "Row Facet Var:", c("None", "treat", varNames, "phase", "phase2", "phase5"))
                           ),
                           
                           column(3,
                                  sliderInput("height", "height", min = 100, max = 2000, value = 500),
                                  sliderInput("width", "width", min = 100, max = 2000, value = 1500),
                                  numericInput("res", "Scaling (click Make Graph after change)", value = 64),
                                  numericInput("dot_size", "Point Size", value = 0.3),
                                  numericInput("gradient_size", "Change density gradient", value = 0.01),
                                  numericInput("text_globe", "Global text size", value = 28),
                                  numericInput("text_axis", "Axis text size", value = 24),
                                  numericInput("text_leg", "Legend text size", value = 20),
                           ),
                           column(2,
                                  radioButtons("filtering", "Filter groups?", c("Yes", "No"), selected = "No"),
                                  checkboxGroupInput("levels", "Groups to include:", c(levels(as.factor(hci.edu.c$treat)))),
                                  radioButtons("phase_filtering", "Filter phases?", c("Yes", "No"), selected = "No"),
                                  checkboxGroupInput("phases", "Phases to include:", c(levels(as.factor(gated.edu.f$phase5))))
                           ),
                           column(3,
                                  selectInput("tran_x", "X-Axis Transformation", c("Linear" = "identity",
                                                                                   "Log10" = "log10",
                                                                                   "pseudolog" = "pseudo_log"),
                                              selected = "identity"),
                                  numericInput("sigma_x", "X-axis Sigma (For pseudolog only)", value = NULL),
                                  selectInput("tran_y", "Y-Axis Transformation", c("Linear" = "identity",
                                                                                   "Log10" = "log10",
                                                                                   "pseudolog" = "pseudo_log"),
                                              selected = "log10"),
                                  numericInput("sigma_y", "Y-axis Sigma (For pseudolog only)", value = NULL),
                                  textInput("title", "Title"),
                                  textInput("x_lab", "X-label"),
                                  textInput("y_lab", "Y-label"),
                                  selectInput("twocolor", "Select Color Map", c("Magma" = "mycolorm",
                                                                                "Rainbow" = "mycolorr",
                                                                                "Viridis" = "mycolorv"),
                                              selected = "mycolorv")),
                           column(2,
                                  numericInput("lower_x", "Lower X Limit", value = NULL),
                                  numericInput("upper_x", "Upper X Limit", value = NULL),
                                  numericInput("lower_y", "Lower Y Limit", value = NULL),
                                  numericInput("upper_y", "Upper Y Limit", value = NULL),
                                  numericInput("aspect", "Aspect Ratio (No Number for free axis)", value = NULL),
                                  numericInput("y_gate", "Y-Axis gate threshold", value = NULL),
                           ))

                           
            ),
    tabPanel("Run Gating", fluidRow(column(6,
                           selectInput("dataset_g", "Select Data Set", c( "Gated" = "gated.edu.f",
                                                                        "Telophases" = "v.unpair",
                                                                        "Late Telophase" = "late_telophase")),
                            numericInput("multi", "Multiplier of Gate", value = 1),
                            selectInput("groups", "Grouping variables (fold_test must go last, don't select phase unless the gated dataframe is selected)", c("treat", varNames, "phase5", "fold_test"), multiple = TRUE),
                            actionButton("fold_test", "Run Gating"),
                            dataTableOutput("DF")))),
    tabPanel("Plot Percents",
             fluidRow(column(12,
                             plotOutput("plot_b", inline = TRUE))),
             fluidRow(column(2,
                             actionButton("graph_click_b", "Make Graph"),
                             varSelectInput("x_var_b", "X_Axis Variable", dplyr::select(hci.edu.c,where(is.factor)),
                                            selected = "cor_int_dnacor"),
                             radioButtons("faceting_b", "Facet Type", c("None", "Wrap", "Grid")),
                             selectInput("colFacet_b", "Column/Wrap Facet Var:", c("None", "treat", varNames, "phase5")),
                             selectInput("colRow_b", "Row Facet Var:", c("None", "treat", varNames, "phase5"))
             ),
             
             column(3,
                    sliderInput("height_b", "height", min = 100, max = 2000, value = 500),
                    sliderInput("width_b", "width", min = 100, max = 2000, value = 1500),
                    numericInput("res_b", "Scaling (click Make Graph after change)", value = 64),

                    numericInput("text_globe_b", "Global text size", value = 28),
                    numericInput("text_axis_b", "Axis text size", value = 24),
                    numericInput("text_leg_b", "Legend text size", value = 20),
                    numericInput("text_bars", "Size of text in bars", value = 6),
                    numericInput("angle", "X-text Angle", value = 45),
                    numericInput("hjust", "X-text alignment", value = 1)
             ),
             column(2,
                    radioButtons("filtering_b", "Filter groups?", c("Yes", "No"), selected = "No"),
                    checkboxGroupInput("levels_b", "Groups to include:", c(levels(as.factor(hci.edu.c$treat)))),
                    radioButtons("phase_filtering_b", "Filter phases?", c("Yes", "No"), selected = "No"),
                    checkboxGroupInput("phases_b", "Phases to include:", c(levels(as.factor(gated.edu.f$phase5))))
             ),
             column(3,

                    textInput("title_b", "Title"),
                    textInput("x_lab_b", "X-label"),
                    textInput("y_lab_b", "Y-label"),
                    textInput("col_lab", "Color-label"),
             column(3,
                    numericInput("aspect_b", "Aspect Ratio (No Number for free axis)", value = NULL)
             ))),
             fluidRow(column(3,
                             textInput("graph_file", "Graph File Name")),
                      column(2,
                             downloadButton("graph", "Save Graph (Use .pdf extensions):")))     
           
        )
      )
    )
    
    ,
    
    server <- function(input, output, session) {
      gg <- reactive(  
        get(input$dataset) %>% 
          {if(input$filtering == "Yes") filter(.,treat %in% c(input$levels)) else .}  %>%
          {if(input$phase_filtering == "Yes") filter(.,phase5 %in% c(input$phases)) else .}  %>% 
          ggplot(aes(x=!!input$x_var, y=!!input$y_var))+
          geom_pointdensity(size = input$dot_size, 
                            adjust = input$gradient_size) +
          scale_y_continuous(transform = if(input$tran_y == "pseudo_log"){scales::pseudo_log_trans(sigma = input$sigma_y)}
                             else{input$tran_y},
                             limits = c(input$lower_y, 
                                        input$upper_y))+
          scale_x_continuous(transform = if(input$tran_x == "pseudo_log"){scales::pseudo_log_trans(sigma = input$sigma_x)}
                             else{input$tran_x},
                             limits = c(input$lower_x, 
                                        input$upper_x))+
          scale_color_gradientn(colours = get(input$twocolor))+
          geom_hline(yintercept = input$y_gate)+
          labs(title = input$title,
               x= input$x_lab,
               y=input$y_lab, 
               colour = "Density")+
          theme_classic()+
          theme(text = element_text(size=input$text_globe),
                legend.title = element_text(size = input$text_leg),
                axis.text = element_text(size=input$text_axis),
                plot.title = element_text(hjust = 0.5),
                axis.line = element_blank(),
                strip.background = element_blank(),
                panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
          {if(is.na(input$aspect)) theme(aspect.ratio=NULL) else theme(aspect.ratio = input$aspect)}+
          {if(input$faceting == "Grid")
            facet_grid(cols = eval(if(input$colFacet == "None") NULL else vars(get(input$colFacet))),
                       rows = eval(if(input$colRow == "None") NULL else vars(get(input$colRow))))}+
          {if(input$faceting == "Wrap")
            facet_wrap(eval(if(input$colFacet == "None") NULL else vars(get(input$colFacet))))} 
      )
      
          loaded <- reactiveValues(data = NULL)
           observeEvent(input$fold_test, {
                        loaded$data <- get(input$dataset_g) %>% foldtest(.,input$y_var, input$y_gate, "fold", input$multi) %>% percents(., groups = input$groups)
                        output$DF <- renderDataTable(loaded$data)})
      
      bar <- reactive(  
        loaded$data %>% 
          {if(input$filtering_b == "Yes") filter(.,treat %in% c(input$levels_b)) else .}  %>%
          {if(input$phase_filtering_b == "Yes") filter(.,phase5 %in% c(input$phases_b)) else .}  %>% 
          ggplot(aes(x=!!input$x_var_b, y=per))+
          geom_col(aes(fill=fct_rev(fold_test)), color = "black", linewidth = 1, width = 0.75)+
          scale_fill_manual(values = c("darkorange1","deepskyblue"))+
          labs(title = input$title_b,
               x= input$x_lab_b,
               y=input$y_lab_b, 
               fill = input$col_lab)+
          theme_classic()+
          theme(text = element_text(size=input$text_globe_b),
                legend.title = element_text(size = input$text_leg_b),
                axis.text.y = element_text(size=input$text_axis_b),
                axis.text.x = element_text(size=input$text_axis_b, angle = input$angle, hjust=input$hjust),
                plot.title = element_text(hjust = 0.5),
                axis.line = element_blank(),
                strip.background = element_blank(),
                panel.background = element_blank())+
          geom_text(
            aes(label = ifelse(round(per,1) >1.9, paste(round(per,1)), paste("")), y = per),
            position = position_stack(vjust = 0.2),
            vjust = 0,
            fontface='bold',
            size=input$text_bars)+
          {if(is.na(input$aspect_b)) theme(aspect.ratio=NULL) else theme(aspect.ratio = input$aspect_b)}+
          {if(input$faceting_b == "Grid")
            facet_grid(cols = eval(if(input$colFacet_b == "None") NULL else vars(get(input$colFacet_b))),
                       rows = eval(if(input$colRow_b == "None") NULL else vars(get(input$colRow_b))))}+
          {if(input$faceting_b == "Wrap")
            facet_wrap(eval(if(input$colFacet_b == "None") NULL else vars(get(input$colFacet_b))))} 
      )
      
      
      observeEvent(input$graph_click,{
        output$plot <- renderPlot({
          gg()
        },
        width = function() input$width,
        height = function() input$height,
        res = input$res
        )
      })
      
      observeEvent(input$graph_click_b,{
        output$plot_b <- renderPlot({
          bar()
        },
        width = function() input$width_b,
        height = function() input$height_b,
        res = input$res_b
        )
      })
      
      
      output$graph <- downloadHandler(filename = function(){
        input$graph_file}, 
        content =  function(file) {
          ggsave(file, gg(), device = "pdf",
                 width = input$width,
                 height = input$height,
                 units = c("px"),
                 dpi = input$res) 
          #pdf(file=file)
          #plot(gg()) 
          #dev.off()
        }
      )
    }
  )
}


