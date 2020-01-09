multiplexR_ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("hr {border-top: 2px solid #000000;}"))
  ),
  
  titlePanel(title="MultiplexR",
             windowTitle="MultiplexR"),
  
  sidebarLayout(
    sidebarPanel("Plot parameters",
      
      fileInput(inputId = "user_file",
                label = "Multiplex data excel or csv file"),
      
      conditionalPanel(
        condition="!!output.",
        selectInput(inputId="cytokine", 
                  label="Cytokine:",
                  choices=c("Cylinders" = "cyl",
                               "Transmission" = "am",
                               "Gears" = "gear")),
                 tableOutput("data"),
    
    selectInput(inputId="bio_rep", 
                label="N3 or N2?:",
                choices=c("Cylinders" = "cyl",
                          "Transmission" = "am",
                          "Gears" = "gear")),
    tableOutput("data"),
    
    selectInput(inputId="statistical_test", 
                label="Statistical Test:",
                choices=c("Student's t-test (unequal variance; unpaired)" = "t-test_unequalVar_unpaired",
                          "Student's t-test (equal variance; unpaired)" = "t-test_equalVar_unpaired",
                          "Student's t-test (unequal variance; paired)" = "t-test_unequalVar_paired",
                          "Student's t-test (equal variance; paired)" = "t-test_equalVar_paired",
                          "Mann Whitney U test" = "Mann_Whitney"),
                selected = "Student's t-test (unequal variance; unpaired)"),
    tableOutput("data"),
    radioButtons(inputId="p_val_display", 
                label="p-value display",
                choices=c("numeric value" = "p_val_num",
                          "stars" = "pval_star"),
                selected = "numeric value"),
    tableOutput("data"),
    
    radioButtons(inputId="show_means", 
                label="Show mean values on plot?",
                choices=c("yes" = "TRUE",
                          "no" = "FALSE"),
                selected = "no"),
    tableOutput("data")
    
    )),
    
    mainPanel("Cytokine comcentration plot",
              plotOutput("Cytokines concentration plot"),
              br(),br(),
              tableOutput("results")
              )
  )
  
)