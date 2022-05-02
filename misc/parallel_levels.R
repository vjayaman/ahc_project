source("result_processing/rfuncs.R")


num_vecs <- 15
file_index <- 3

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
p <- 4
accuracy_tbl <- parallelContainedElementsTable(s_x3, p) %>% 
  set_colnames(c("level", paste0("par_0_elements")))

proc_p <- containedElementsTable(s_x3[proc == p]) %>% 
  set_colnames(c(paste0("proc_",p,"_level"), paste0("par_", p, "_elements")))

# now just need to add the elements from proc 4 to any element list of 
# accuracy table that contains the parent of those children







