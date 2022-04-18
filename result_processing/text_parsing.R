# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/v1/")

library(magrittr)
library(tidyverse)
# library(ggplot2)
library(fastcluster)

# # PARALLEL ---------------------------------------------------------------------
# 
# pfiles <- list.files("outputs/parallel/", full.names = TRUE)
# pnames <- pfiles %>% gsub("outputs/parallel/", "", .) %>% gsub(".txt", "", .)
# names(pfiles) <- pnames
# 
# 
# parallel_times <- lapply(pfiles, function(x) {
#   p_outputs <- readLines(x)
#   p_ind1 <- which(grepl("Total time", p_outputs))
#   y <- strsplit(x, "/")[[1]][3] %>% gsub(".txt", "", .) %>% 
#     strsplit(names(x), split = "m") %>% unlist()
#   z <- strsplit(p_outputs[p_ind1], " ")[[1]][4] %>% as.double()
#   data.frame(t(y), z) %>% set_colnames(c("num_nodes", "num_vecs", "time"))
# }) %>% bind_rows()
# 
# 
# # NAIVE SEQUENTIAL -------------------------------------------------------------
# 
sfiles <- list.files("outputs/seq/", full.names = TRUE)
# 
# naive_sequential_times <- lapply(sfiles, function(x) {
#   s_outputs <- readLines(x)
#   s_ind1 <- which(grepl("Total time", s_outputs))
#   y <- strsplit(x, "/")[[1]][3] %>% gsub(".txt", "", .) %>% 
#     strsplit(names(x), split = "m") %>% unlist() %>% extract2(2)
#   z <- strsplit(s_outputs[s_ind1], " ")[[1]][4] %>% as.double()
#   data.frame(num_nodes = 0, num_vecs = as.double(y), time = z)
# }) %>% bind_rows()
# 

# SEQUENTIAL -------------------------------------------------------------------

num_vec_list <- lapply(sfiles, function(sf) {
  strsplit(sf, "/")[[1]][3] %>% gsub(".txt", "", .) %>% gsub("m", "", .) %>% as.double()
}) %>% unlist()


sequential_times <- lapply(1:length(num_vec_list), function(i) {
  file_index <- i
  num_vecs <- num_vec_list[i]
  
  stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)
  
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
  hc <- fastcluster::hclust(df, members = "ave")
  # plot(hc)
  # plot(hc, hang = -1)
  
  stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
  
  df <- data.frame(num_nodes = 0, num_vecs, 
                   time = difftime(stopwatch[['end_time']], stopwatch[['start_time']], units = "secs") %>% as.double())
  # print(df)
  # df
}) %>% bind_rows()



