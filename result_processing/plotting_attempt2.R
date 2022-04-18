source("rfuncs.R")


num_vecs <- 5 #15
file_index <- 1 #2

# ------------------------------------------------------------------------------
# SEQUENTIAL -------------------------------------------------------------------

# # NOW FOR PLOTTING
# 
# base_s <- s_x3 %>% arrange(-level) %>% select(lvl, seqs) %>% unique() %>% 
#   rownames_to_column(var = "x") %>% add_column(y = 0) %>% mutate(across(x, as.double))
# all_levels <- base_s[lvl <= 2]
# p <- ggplot() +
#   geom_point(data = all_levels[lvl == 1], aes(x = x, y = y)) + ylim(0, nrow(base_s))
# 
# go_through <- s_x3[variable == "center"]
# 
# for (lvl_x in go_through$lvl) {
#   ax <- go_through[lvl == lvl_x]
#   for (m in 1:nrow(ax)) {
#     a1 <- ax[m,]
#     b1 <- s_x3[level == a1$level][variable != "center"] %>% 
#       arrange(seqs) %>%  filter(row_number()==1 | row_number()==n())
#     prev_lvl <- b1[lvl != lvl_x] %>% select(seqs, lvl)
#     c1 <- all_levels %>% inner_join(prev_lvl, by = c("seqs", "lvl"))
#     
#     if (nrow(c1) == 0) {
#       c1 <- base_s %>% inner_join(prev_lvl, by = c("seqs", "lvl"))
#     }else if (nrow(c1) == 1) {
#       d1 <- base_s[seqs == prev_lvl$seqs[!(prev_lvl$seqs %in% c1)] & 
#                      lvl == prev_lvl$lvl[!(prev_lvl$seqs %in% c1)]]
#       c1 <- bind_rows(c1, d1)
#     }else if (nrow(c1) > 2) {
#       print(paste0("problem found for lvl ", lvl_x, ", m = ", m))
#     }
#     c1$y <- c1$lvl
#     
#     if (nrow(c1) == 2) {
#       a1$x <- sum(c1$x)/2
#       a1$y <- max(c1$y) + 1
#       all_levels <- a1 %>% select(x, lvl, seqs, y) %>% bind_rows(all_levels, .)
#       p <- p + geom_point(data = a1, aes(x = x, y = y))  
#     }
#   }
# }



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


p_x2 <- lapply(p_pairs, function(x_i) {
  x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist() %>% t() %>% as.data.frame() %>% 
    select(-1)
}) %>% bind_rows() %>% set_colnames(c("x1", "x2", "center")) %>% 
  add_column(level = p_levels$level, .before = 1) %>% 
  add_column(proc = p_levels$proc, .before = 1) %>% 
  sepCols(3, .) %>% sepCols(4, .) %>% sepCols(5, .) %>% as_tibble()


p_x1 <- lapply(p_pairs, function(x_i) {
  x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
})


p_first_pair <- orderedPairFormatting(p_x1, 1) %>% add_column(type = "1")
p_second_pair <- orderedPairFormatting(p_x1, 2) %>% add_column(type = "1")
p_new_center <- orderedPairFormatting(p_x1, 3) %>% add_column(type = "2")

# # start of plotting ------------------------------------------------------------
# p_toplot <- bind_rows(p_first_pair, p_second_pair) %>% bind_rows(., p_new_center)
# 
# p_d1 <- p_toplot[p_toplot$type == 1,] %>% arrange(val) %>% 
#   add_column(new_val = rep(c(0,1), nrow(p_toplot[p_toplot$type == 1,])/2))
# 
# p_d2a <- p_d1[p_d1$new_val == 0,c("val","x","y","type")]
# p_d2b <- p_d1[p_d1$new_val == 1,c("val","x","y","type")] %>% 
#   set_colnames(c("val","xend","yend","type"))
# p_d2 <- inner_join(p_d2a, p_d2b, by = c("val","type"))
# 
# p_df <- bind_rows(p_d2, p_d2) %>% bind_rows(p_d2)
# 
# p_toplot$type[p_toplot$type == "1"] <- "Merged pair"
# p_toplot$type[p_toplot$type == "2"] <- "Center"
# p_plot <- 
#   ggplot(data = p_toplot, aes(x = x, y = y, color = type)) + geom_point() + 
#   geom_segment(aes(x = p_df$x, y = p_df$y, xend = p_df$xend, yend = p_df$yend)) + 
#   ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
#   labs(color = "Point type")
# p_plot
# ggsave(paste0("parallel_plot.png"))
# # end of plotting --------------------------------------------------------------








