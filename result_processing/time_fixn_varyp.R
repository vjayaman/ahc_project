# setwd("C:/Users/vasen/Documents/MSc courses/7850-APC/ProjectNotes/ahc_project")

library(magrittr)
library(tidyverse)
library(fastcluster)
library(data.table)
library(xtable)

# # PARALLEL: fix n vary p -----------------------------------------------------

pfiles <- list.files("outputs/fix_n_vary_p/parallel/", full.names = TRUE)
pnames <- pfiles %>% gsub("outputs/fix_n_vary_p/parallel/", "", .) %>% gsub(".txt", "", .)
names(pfiles) <- pnames

parallelFixNVaryP <- function(pfiles, m_val) {
  pfiles <- grep(m_val, pfiles, value = TRUE)
  parallel_times <- lapply(1:length(pfiles), function(i) {
    x <- pfiles[i]
    p_outputs <- readLines(x)
    p_ind1 <- which(grepl("Total time", p_outputs))
    if (length(p_ind1) > 0) {
      y <- strsplit(names(x), split = "m") %>% unlist() %>% gsub("p", "", .) %>% as.double()
      z <- strsplit(p_outputs[p_ind1], " ")[[1]][4] %>% as.double()
      data.frame(t(y), z) %>% set_colnames(c("num_nodes", "num_vecs", "time")) 
    }
  }) %>% bind_rows()
  return(parallel_times)
}

# fix n = 15 vary p
# fix n = 60 vary p
# fix n = 300 vary p
parallel_times <- parallelFixNVaryP(pfiles, "m015") %>% 
  bind_rows(., parallelFixNVaryP(pfiles, "m060")) %>% 
  bind_rows(., parallelFixNVaryP(pfiles, "m300")) %>% as.data.table()


# # General seq. -------------------------------------------------------------

sfiles <- list.files("outputs/fix_n_vary_p/seq/", full.names = TRUE)
snames <- sfiles %>% gsub("outputs/fix_n_vary_p/seq/", "", .) %>% gsub(".txt", "", .)
names(sfiles) <- snames

seqFixNVaryP <- function(sfiles, m_val) {
  sfiles <- grep(m_val, sfiles, value = TRUE)
  naive_sequential_times <- lapply(1:length(sfiles), function(i) {
    x <- sfiles[i]
    s_outputs <- readLines(x)
    s_ind1 <- which(grepl("Total time", s_outputs))
    y <- strsplit(names(x), split = "m") %>% unlist() %>% gsub("p", "", .) %>% as.double()
    z <- strsplit(s_outputs[s_ind1], " ")[[1]][4] %>% as.double()
    data.frame(num_nodes = y[1], num_vecs = y[2], time = z)
  }) %>% bind_rows()
  return(naive_sequential_times)
}

naive_sequential_times <- seqFixNVaryP(sfiles, "m015") %>% 
  bind_rows(., seqFixNVaryP(sfiles, "m060")) %>% 
  bind_rows(., seqFixNVaryP(sfiles, "m300")) %>% as.data.table()


# SEQUENTIAL -------------------------------------------------------------------


num_vec_list <- c(15, 60, 300)

sequential_times <- lapply(1:length(num_vec_list), function(i) {
  file_index <- i
  num_vecs <- num_vec_list[i]
  
  start_time <- Sys.time()
  
  sfiles <- list.files("outputs/fix_n_vary_p/seq/", full.names = TRUE)
  snames <- sfiles %>% gsub("outputs/fix_n_vary_p/seq/", "", .) %>% gsub(".txt", "", .)
  names(sfiles) <- snames
  s_outputs <- readLines(sfiles[file_index])
  
  s_ind1 <- which(grepl("Matrix A: ", s_outputs))+1
  seq_results <- s_outputs[s_ind1:(s_ind1+(num_vecs-1))]
  
  s_a1 <- lapply(seq_results, function(x) {
    strsplit(x, "\t") %>% unlist() %>% as.double()
  }) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
    as_tibble() %>% set_colnames(c("x", "y"))
  
  df <- stats::dist(s_a1)
  hc <- fastcluster::hclust(df, members = "ave") # plot(hc); plot(hc, hang = -1)
  
  end_time <- Sys.time() - start_time
  
  data.frame(num_nodes = 1, num_vecs, time = as.character(end_time) %>% as.double())
}) %>% bind_rows() %>% as.data.table()


# RUNNING ----------------------------------------------------------------------
# fix n = 15, vary p
nv <- 15

# a1 <- sequential_times[num_vecs == nv] %>% add_column(type = "Best seq.")
a1 <- sequential_times[num_vecs == nv] %>% 
  set_colnames(c("seq p", "# vectors", "Best seq."))

a2 <- naive_sequential_times[num_vecs == nv] %>% 
  set_colnames(c("seq p", "# vectors", "General seq. (mine)"))

a3 <- parallel_times[num_vecs == nv] %>% 
  set_colnames(c("parallel p", "# vectors", "Parallel (mine)"))

n15 <- merge.data.table(a1, a2, by = c("seq p", "# vectors")) %>% 
  merge.data.table(., a3, by = "# vectors") %>% 
  mutate(Speedup = `Best seq.` / `Parallel (mine)`)

# n15 <- bind_rows(a1, a2) %>% bind_rows(., a3) %>% 
#   set_colnames(c("Number of processes", "Number of vectors", "Time (sec)", "Implementation"))
print(xtable(x = n15, type = "latex", digits = c(0,0,0,5,5,0,5,5)), include.rownames = FALSE)

# fix n = 60, vary p
nv <- 60
a1 <- sequential_times[num_vecs == nv] %>% 
  set_colnames(c("seq p", "# vectors", "Best seq."))

a2 <- naive_sequential_times[num_vecs == nv] %>% 
  set_colnames(c("seq p", "# vectors", "General seq. (mine)"))

a3 <- parallel_times[num_vecs == nv] %>% 
  set_colnames(c("parallel p", "# vectors", "Parallel (mine)"))

n60 <- merge.data.table(a1, a2, by = c("seq p", "# vectors")) %>% 
  merge.data.table(., a3, by = "# vectors") %>% 
  mutate(Speedup = `Best seq.` / `Parallel (mine)`)

print(xtable(x = n60, type = "latex", digits = c(0,0,0,5,5,0,5,5)), include.rownames = FALSE)

# fix n = 300, vary p
nv <- 300

a1 <- sequential_times[num_vecs == nv] %>% 
  set_colnames(c("seq p", "# vectors", "Best seq."))

a2 <- naive_sequential_times[num_vecs == nv] %>% 
  set_colnames(c("seq p", "# vectors", "General seq. (mine)"))

a3 <- parallel_times[num_vecs == nv] %>% 
  set_colnames(c("parallel p", "# vectors", "Parallel (mine)"))

n300 <- merge.data.table(a1, a2, by = c("seq p", "# vectors")) %>% 
  merge.data.table(., a3, by = "# vectors") %>% 
  mutate(Speedup = `Best seq.` / `Parallel (mine)`)

print(xtable(x = n300, type = "latex", digits = c(0,0,0,5,5,0,5,5)), include.rownames = FALSE)

