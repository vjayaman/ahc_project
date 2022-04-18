source("result_processing/rfuncs.R")
# ------------------------------------------------------------------------------
# SEQUENTIAL -------------------------------------------------------------------

sfiles <- list.files("outputs/seq/", full.names = TRUE)
snames <- sfiles %>% gsub("outputs/seq/", "", .) %>% gsub(".txt", "", .)
names(sfiles) <- snames

num_vecs <- 15
file_index <- 3

s_outputs <- readLines(sfiles[file_index])

s_ind1 <- which(grepl("Matrix A: ", s_outputs))+1
seq_results <- s_outputs[s_ind1:(s_ind1+(num_vecs-1))]

s_a1 <- initialPairs(seq_results)

s_lvlinds <- which(grepl("Minimum distance ", s_outputs))
s_x2 <- levelCharNumTbl(s_lvlinds, s_outputs)

# plotting
 
s_x3 <- buildLvlTbl(s_x2, s_a1, FALSE)$s_x3

accuracy_tbl <- containedElementsTable(s_x3) %>% set_colnames(c("level", "seq_elements"))

