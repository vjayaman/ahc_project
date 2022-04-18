source("result_processing/rfuncs.R")


num_vecs <- 10
file_index <- 2

# ------------------------------------------------------------------------------
# PARALLEL ---------------------------------------------------------------------

pfiles <- list.files("outputs/parallel/", full.names = TRUE)
pnames <- pfiles %>% gsub("outputs/parallel/", "", .) %>% gsub(".txt", "", .)
names(pfiles) <- pnames
p_outputs <- readLines(pfiles[file_index])

p_ind1 <- which(grepl("Matrix A: ", p_outputs))[1]+1
p_matA <- p_outputs[p_ind1:(p_ind1+(num_vecs-1))]

p_a1 <- lapply(p_matA, function(x) {
  strsplit(x, "\t") %>% unlist() %>% as.double()
}) %>% as.data.frame() %>% t() %>% as.data.frame() %>% 
  as_tibble() %>% set_colnames(c("x", "y"))


# p_lvlinds <- which(grepl("Minimum distance ", p_outputs))
p_pairs <- grep("From ", p_outputs, value = TRUE)
p_clustering <- grep("Minimum distance ", p_outputs, value = TRUE)

# note that the "levels" are decreasing --> max is at the leaves, 1 is at the root

p_levels <- lapply(p_clustering, function(x_i) {
  y_i <- strsplit(strsplit(x_i, split = "level ")[[1]][2], ", process ")[[1]] %>% as.double()
  data.table(level = y_i[1], proc = y_i[2])
}) %>% bind_rows()


df <- lapply(p_pairs, function(x_i) {
  y_i <- x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
  y_i[nchar(y_i) > 0] %>% t() %>% as.data.frame()
}) %>% bind_rows() 

# # for m = 15
# df[7,1] <- df[7,4]
# df[7,2] <- df[7,5]
# df[7,3] <- df[7,6]
# df <- df[,1:3]

p_x2 <- df %>% set_colnames(c("x1", "x2", "center")) %>% 
  add_column(level = p_levels$level, .before = 1) %>% 
  add_column(proc = p_levels$proc, .before = 1) %>% 
  sepCols(3, .) %>% sepCols(4, .) %>% sepCols(5, .) %>% as_tibble() %>% unique() %>% 
  as.data.table()


p_x1 <- lapply(p_pairs, function(x_i) {
  x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
})

p_x2[proc == 0]$level <- p_x2[proc == 0]$level + max(p_x2[proc != 0]$level)
s_x3 <- buildLvlTbl(p_x2, p_a1, TRUE)$s_x3


# need to run a version of containedElementsTable() that considers proc 0 
# and one other process at a time. 
# i.e. when tracing other sequences, starting from proc 0, 
# if they are found as an internal node for a proc that is not the given 
# proc, then consider it as a leaf and stop the recursion
# then, should return a separate table for each worker process (that includes)
# the clustering done closer to the root at proc 0

# for (p in sort(unique(p_x2$proc))[-1]) {
#   # s_x3 <- s_x3[proc %in% c(0,1)]
#   accuracy_tbl <- containedElementsTable(s_x3) %>% 
#     set_colnames(c("level", "par_proc0_elements"))
#   print(paste0("Process ", p))
#   print(accuracy_tbl)
# }
parallelContainedElementsTable(s_x3, 4)


# generate table of level and absorbed elements into a group
# designed in a sequential way
# will use this for comparing "accuracy" for sequential and parallel groups later
parallelContainedElementsTable <- function(s_x3, p) {
  dif_levels <- sort(unique(s_x3$level), decreasing = TRUE)
  for (r in dif_levels) {
    # print(r)
    if (r == max(unique(s_x3$level))) {
      accuracy_tbl <- data.table(
        level = r, elements = list(s_x3[level == r & variable != "center"]$value))
    }else {
      merged_pair <- s_x3[level == r & variable != "center"]
      leaves <- merged_pair[lvl == 1]
      if (nrow(leaves) > 0) {
        seq_list <- list(leaves$value)
      }else {
        seq_list <- list()
      }
      # a pair merged at this level that have children
      x1 <- merged_pair[lvl != 1]
      if (nrow(x1) > 0) {
        for (i in 1:nrow(x1)) {
          seq_list <- append(seq_list, parTraceSequences(s_x3, x1[i,], r, seq_list, p))
          df <- tibble(level = r, elements = list(unique(unlist(seq_list))))
          accuracy_tbl <- bind_rows(accuracy_tbl, df)
        }  
      }else { # just merging two leaves
        accuracy_tbl <- bind_rows(
          accuracy_tbl, 
          data.table(level = r, elements = list(s_x3[level == r & variable != "center"]$value)))
      }
    }
  }
  return(accuracy_tbl)
}


parTraceSequences <- function(s_x3, x, r, seq_list, p) {
  if (nrow(x) > 0) {
    for (j in 1:nrow(x)) {
      z <- x[j,]
      # when was this value formed in a merge?
      # i.e. when it was a center, at a level before the cuurrent one, 
      # and when the sequence ID was the same
      # level_of_earlier_merge <- s_x3[variable == "center" & lvl < z$lvl & seqs == z$seqs]$level
      
      # arbitrarily pick the first of these (since we consider it a leaf if proc != given or 0)
      earlier_merge <- s_x3[seqs == z$seqs & variable == "center"][1]
      
      # collect the triple of x1, x2, pair at this earlier level
      y <- s_x3[level == earlier_merge$level & proc == earlier_merge$proc]
      
      # what was the merged pair?
      mp <- y[variable != "center"]
      
      # add any children from that pair to the list if they were leaves or from a different proc
      seq_list <- append(seq_list, mp[lvl == 1 & proc %in% c(0, p)]$value)
      seq_list <- append(seq_list, mp[lvl == 1 & !(proc %in% c(0,p))]$value)
      # if they were not leaves, then we want to trace back to find when they 
      # were merged into an internal node -- will keep doint this until we have 
      # added all leaves to the list
      z <- mp[lvl != 1 & proc %in% c(0, p)]
      seq_list <- append(seq_list, parTraceSequences(s_x3, z, r, seq_list, p))
    }
  }else {
    # print(list(unlist(seq_list)))
    return(list(unlist(seq_list)))
  }
  return(seq_list)
}








