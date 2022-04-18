source("rfuncs.R")
# ------------------------------------------------------------------------------
# SEQUENTIAL -------------------------------------------------------------------

# collectSeqs()
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


s_lvlinds <- which(grepl("Minimum distance ", s_outputs))
s_clustering <- s_outputs[s_lvlinds+1]

s_levels <- lapply(s_outputs[s_lvlinds], function(x_i) {
  strsplit(x_i, split = "level ")[[1]][2] %>% as.double()
}) %>% unlist()

s_x2 <- lapply(s_clustering, function(x_i) {
  x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist() %>% t() %>% as.data.frame() %>% 
    select(-1)
}) %>% bind_rows() %>% set_colnames(c("x1", "x2", "center")) %>% 
  add_column(level = s_levels, .before = 1) %>% 
  sepCols(2, .) %>% sepCols(3, .) %>% sepCols(4, .) %>% as_tibble()

# plotting
 
s_x3 <- buildLvlTbl(s_x2, s_a1)$s_x3
for (r in unique(s_x3$level)) {
  if (r == num_vecs) {
    accuracy_tbl <- data.table(
      level = r, elements = list(s_x3[level == r & variable != "center"]$value))
  }else {
    merged_pair <- s_x3[level == r & variable != "center"]
    seq_list <- list(merged_pair[lvl == 1]$value)
    x1 <- merged_pair[lvl != 1]
    for (i in 1:nrow(x1)) {
      df <- traceSequences(s_x3, x1[i,], r, seq_list)
      accuracy_tbl <- bind_rows(accuracy_tbl, df)
    }  
  }
}

accuracy_tbl <- accuracy_tbl %>% set_colnames(c("level", "seq_elements"))
# accuracy_tbl

