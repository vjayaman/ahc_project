library(shiny)
library(shinydashboard)
library(DT)
source("result_processing/rfuncs.R")

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    selectInput(inputId = "num_vecs",
      "Number of vectors",
      c("0005", "0010", "0015", "0025", "0050", "0075", "0100", "0250", "0500", "0750", "1000", "1500", "2000", "3000"),
      selected = 1, multiple = FALSE), 
    
    sidebarMenu(
      menuItem(text = "Fix p and vary n", tabName = "fix_p"), 
      # menuItem(text = "Fix n and vary p", tabName = "fix_n"), # did this manually for report
      menuItem(text = "General computation times", tabName = "times"), 
      menuItem(text = "Efficiency", tabName = "eff")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem("fix_p", 
              fluidRow(
                box(plotOutput("parallel_plot_halfway", height = 350)), 
                box(plotOutput("parallel_plot", height = 350))
              ), 
              fluidRow(
                box(plotOutput("seq_plot_halfway", height = 350)), 
                box(plotOutput("seq_plot", height = 350))
              )
            ),
      # tabItem("fix_n", "Fix n plot comparison"), # did this manually for report
      tabItem("times", 
              fluidRow(box(title = "Computation times for parallel implementation", 
                           width = 4, DT::dataTableOutput("parallel_times")), 
                       box(title = "Computation times for sequential implementation", 
                           width = 4, DT::dataTableOutput("naive_seq_times")), 
                       box(title = "Computation times for running fastcluster (Mullner)", 
                           width = 4, DT::dataTableOutput("best_seq_times"))
              ), 
              fluidRow(box(title = "Computation times for all versions, compared", 
                           DT::dataTableOutput("all_times")))), 
      tabItem("eff", 
              fluidRow(box(title = "Fixed number of vectors at 15, varied p", DT::dataTableOutput("fix_n15"))), 
              fluidRow(box(title = "Fixed number of vectors at 60, varied p", DT::dataTableOutput("fix_n60"))), 
              fluidRow(box(title = "Fixed number of vectors at 300, varied p", DT::dataTableOutput("fix_n300"))))
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
  
  dfs <- reactiveValues()
  times <- reactiveValues(nst = FALSE, st = FALSE, pt = FALSE)
  
  # Fix p and vary n
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
  
  # Fix n and vary p
  # General computation times
  output$parallel_times <- DT::renderDataTable({
    parallel_times <- lapply(pfiles, function(x) {
      p_outputs <- readLines(x)
      p_ind1 <- which(grepl("Total time", p_outputs))
      if (length(p_ind1) > 0) {
        y <- strsplit(x, "/")[[1]][3] %>% gsub(".txt", "", .) %>%
          strsplit(names(x), split = "m") %>% unlist() %>% gsub("p", "", .) %>% 
          as.double()
        z <- strsplit(p_outputs[p_ind1], " ")[[1]][4] %>% as.double()
        data.frame(t(y), z) %>% set_colnames(c("num_nodes", "num_vecs", "time")) 
      }
    }) %>% bind_rows() %>% 
      set_colnames(c("# processes", "# vectors", "Time (sec)"))
    dfs$par <- parallel_times
    
    DT::datatable(parallel_times, rownames = FALSE)
  })
  
  output$naive_seq_times <- DT::renderDataTable({
    naive_sequential_times <- lapply(sfiles, function(x) {
      s_outputs <- readLines(x)
      s_ind1 <- which(grepl("Total time", s_outputs))
      y <- strsplit(x, "/")[[1]][3] %>% gsub(".txt", "", .) %>%
        strsplit(names(x), split = "m") %>% unlist() %>% extract2(2)
      z <- strsplit(s_outputs[s_ind1], " ")[[1]][4] %>% as.double()
      data.frame(num_nodes = 1, num_vecs = as.double(y), time = z)
    }) %>% bind_rows() %>% 
      set_colnames(c("# processes", "# vectors", "Time (sec)"))
    dfs$nst <- naive_sequential_times
    
    DT::datatable(naive_sequential_times, rownames = FALSE)
  })

  output$best_seq_times <- DT::renderDataTable({
    num_vec_list <- lapply(sfiles, function(sf) {
      strsplit(sf, "/")[[1]][3] %>% gsub(".txt", "", .) %>% gsub("m", "", .) %>% as.double()
    }) %>% unlist()

    sequential_times <- lapply(1:length(num_vec_list), function(i) {
      file_index <- i
      num_vecs <- num_vec_list[i]
      
      start_time <- Sys.time()
      
      sfiles <- list.files("outputs/seq/", full.names = TRUE)
      snames <- sfiles %>% gsub("outputs/seq/", "", .) %>% gsub(".txt", "", .)
      names(sfiles) <- snames
      s_outputs <- readLines(sfiles[file_index])
      
      s_ind1 <- which(grepl("Matrix A: ", s_outputs))+1
      seq_results <- s_outputs[s_ind1:(s_ind1+(num_vecs-1))]
      
      s_a1 <- lapply(seq_results, function(x) {
        strsplit(x, "\t") %>% unlist() %>% as.double()
      }) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
        as_tibble() %>% set_colnames(c("x", "y"))
      
      df <- stats::dist(s_a1)
      hc <- fastcluster::hclust(df, members = "ave") # plot(hc) # plot(hc, hang = -1)
      
      end_time <- Sys.time() - start_time
      
      df <- data.frame(num_nodes = 1, num_vecs, 
                       time = as.character(end_time) %>% as.double())
    }) %>% bind_rows() %>% 
      set_colnames(c("# processes", "# vectors", "Time (sec)"))
    dfs$st <- sequential_times
    
    DT::datatable(sequential_times, rownames = FALSE)
  })
  
  output$all_times <- DT::renderDataTable({
    pt <- dfs$par %>% select(2,3) %>% set_colnames(c("# vectors", "Parallel (sec)"))
    nst <- dfs$nst %>% select(2,3) %>% set_colnames(c("# vectors", "Naive seq. (sec)"))
    st <- dfs$st %>% select(2,3) %>% set_colnames(c("# vectors", "Best seq. (sec)"))
    
    times <- merge.data.table(pt, nst, by = "# vectors") %>% 
      merge.data.table(., st, by = "# vectors") %>% 
      mutate(Speedup = `Best seq. (sec)` / `Parallel (sec)`, 
             `Cost Optimality` = `Best seq. (sec)` / (as.double(unique(dfs$par$`# processes`)) * `Parallel (sec)`))
    DT::datatable(times, rownames = FALSE)
  })
  
  
  # Efficiency
  output$fix_n15 <- DT::renderDataTable({
    nv <- 15
    
    times$st <- runBestSeq("outputs/fix_n_vary_p/seq/")
    times$nst <- runseqFixNVaryP("outputs/fix_n_vary_p/seq/")
    times$pt <- runParallelFixNVaryP("outputs/fix_n_vary_p/parallel/")
    
    a1 <- times$st[num_vecs == nv] %>% set_colnames(c("seq p", "# vectors", "Best seq."))
    a2 <- times$nst[num_vecs == nv] %>% set_colnames(c("seq p", "# vectors", "General seq. (mine)"))
    a3 <- times$pt[num_vecs == nv] %>% set_colnames(c("parallel p", "# vectors", "Parallel (mine)"))
    
    n15 <- merge.data.table(a1, a2, by = c("seq p", "# vectors")) %>% 
      merge.data.table(., a3, by = "# vectors") %>% 
      mutate(Speedup = `Best seq.` / `Parallel (mine)`)
    
    DT::datatable(n15, rownames = FALSE)
  })
  
  output$fix_n60 <- DT::renderDataTable({
    req(times$st, times$nst, times$pt)
    nv <- 60
    a1 <- times$st[num_vecs == nv] %>% 
      set_colnames(c("seq p", "# vectors", "Best seq."))
    a2 <- times$nst[num_vecs == nv] %>% 
      set_colnames(c("seq p", "# vectors", "General seq. (mine)"))
    a3 <- times$pt[num_vecs == nv] %>% 
      set_colnames(c("parallel p", "# vectors", "Parallel (mine)"))
    
    n60 <- merge.data.table(a1, a2, by = c("seq p", "# vectors")) %>% 
      merge.data.table(., a3, by = "# vectors") %>% 
      mutate(Speedup = `Best seq.` / `Parallel (mine)`)
    
    DT::datatable(n60, rownames = FALSE)
  })
  
  output$fix_n300 <- DT::renderDataTable({
    req(times$st, times$nst, times$pt)
    nv <- 300
    
    a1 <- times$st[num_vecs == nv] %>% 
      set_colnames(c("seq p", "# vectors", "Best seq."))
    a2 <- times$nst[num_vecs == nv] %>% 
      set_colnames(c("seq p", "# vectors", "General seq. (mine)"))
    a3 <- times$pt[num_vecs == nv] %>% 
      set_colnames(c("parallel p", "# vectors", "Parallel (mine)"))
    
    n300 <- merge.data.table(a1, a2, by = c("seq p", "# vectors")) %>% 
      merge.data.table(., a3, by = "# vectors") %>% 
      mutate(Speedup = `Best seq.` / `Parallel (mine)`)
    
    DT::datatable(n300, rownames = FALSE)
  })
  
  
  
}

shinyApp(ui, server)