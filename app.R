library(shiny)
library(shinydashboard)
source("result_processing/rfuncs.R")

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    selectInput(inputId = "num_vecs",
      "Number of vectors",
      c("0005", "0010", "0015", "0025", "0050", "0075", "0100", "0250", "0500", "0750", "1000", "1500", "2000", "3000"),
      selected = 1, multiple = FALSE)
  ),
  dashboardBody(
    fluidRow(
      box(plotOutput("parallel_plot_halfway", height = 350)), 
      box(plotOutput("parallel_plot", height = 350))
    ), 
    fluidRow(
      box(plotOutput("seq_plot_halfway", height = 350)), 
      box(plotOutput("seq_plot", height = 350))
    )
  )
)

server <- function(input, output) {
  pfiles <- list.files("outputs/parallel/", full.names = TRUE)
  pnames <- pfiles %>% gsub("outputs/parallel/", "", .) %>% gsub(".txt", "", .)
  names(pfiles) <- pnames
  sfiles <- list.files("outputs/seq/", full.names = TRUE)
  snames <- sfiles %>% gsub("outputs/seq/", "", .) %>% gsub(".txt", "", .)
  names(sfiles) <- snames
  
  
  output$parallel_plot_halfway <- renderPlot({
    p_outputs <- readLines(pfiles[grep(input$num_vecs, pnames)])
    parPlotsGeneral(as.double(input$num_vecs), p_outputs, plot_level = 1)
  })
  
  output$parallel_plot <- renderPlot({
    p_outputs <- readLines(pfiles[grep(input$num_vecs, pnames)])
    parPlotsGeneral(as.double(input$num_vecs), p_outputs, plot_level = 2)
  })
  
  
  output$seq_plot_halfway <- renderPlot({
    s_outputs <- readLines(sfiles[grep(input$num_vecs, pnames)])
    seqPlotsGeneral(as.double(input$num_vecs), s_outputs, plot_level = 1)
  })
  
  output$seq_plot <- renderPlot({
    s_outputs <- readLines(sfiles[grep(input$num_vecs, pnames)])
    seqPlotsGeneral(as.double(input$num_vecs), s_outputs, plot_level = 2)
  })
  
}

shinyApp(ui, server)