# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/ahc_project/")

source("result_processing/rfuncs.R")

num_vecs <- 15
file_index <- 3


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
  t_i <- strsplit(x_i, split = "level ")[[1]][2]
  u_i <- strsplit(t_i, ", process ")[[1]][1] %>% as.double()
  v_i <- strsplit(strsplit(t_i, ", process ")[[1]][2], "\t")[[1]][1] %>% as.double()
  data.table(level = u_i, proc = v_i)
}) %>% bind_rows()

# for m = 10
# p_pairs[4] <- gsub("3.00\t3.00", "", p_pairs[4])

# for m = 10
# p_pairs[14] <- gsub("5.00\t3.50", "", p_pairs[14])

df <- lapply(p_pairs, function(x_i) {
  y_i <- strsplit(x_i, "From ")[[1]][2] %>% 
    # gsub("From ", "", .) %>% 
    gsub(" and ", "\t", .) %>% 
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
  strsplit(x_i, "From ")[[1]][2] %>% 
    # gsub("From ", "", .) %>% 
    gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
})

p_first_pair <- orderedPairFormatting(p_x1, 1) %>% add_column(type = "1")
p_second_pair <- orderedPairFormatting(p_x1, 2) %>% add_column(type = "1")
p_new_center <- orderedPairFormatting(p_x1, 3) %>% add_column(type = "2")




p_toplot <- bind_rows(p_first_pair, p_second_pair) %>% bind_rows(., p_new_center) %>% 
  as.data.table()
# p_toplot <- p_toplot[!(val == 4)]
# p_toplot <- p_toplot[!(x == 0 & y == 0)]

p_d1 <- p_toplot[p_toplot$type == 1,] %>% arrange(val) %>% 
  add_column(new_val = rep(c(0,1), nrow(p_toplot[p_toplot$type == 1,])/2))

p_d2a <- p_d1[p_d1$new_val == 0,c("val","x","y","type")]
p_d2b <- p_d1[p_d1$new_val == 1,c("val","x","y","type")] %>% 
  set_colnames(c("val","xend","yend","type"))
p_d2 <- inner_join(p_d2a, p_d2b, by = c("val","type"))

p_df <- bind_rows(p_d2, p_d2) %>% bind_rows(p_d2)

p_toplot$type[p_toplot$type == "1"] <- "Merged pair"
p_toplot$type[p_toplot$type == "2"] <- "Center"

# p_plot <- 
#   ggplot(data = p_toplot, aes(x = x, y = y, color = type)) + geom_point() + 
#   geom_segment(aes(x = p_df$x, y = p_df$y, xend = p_df$xend, yend = p_df$yend)) + 
#   ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
#   labs(color = "Point type")
# p_plot
# ggsave(paste0("parallel_plot.png"))

p_toplot <- p_toplot %>% mutate(across(val, as.double))
p_df <- p_df %>% mutate(across(val, as.double))

if (num_vecs == 5) {
  toplot_1 <- as.data.table(p_toplot)[val <= 4]
  df_1 <- as.data.table(p_df)[val <= 4]
  
  toplot_2 <- as.data.table(p_toplot)[val > 4]
  df_2 <- as.data.table(p_df)
}else {
  toplot_1 <- as.data.table(p_toplot)[val <= floor(num_vecs/2)]
  df_1 <- as.data.table(p_df)[val <= floor(num_vecs/2)]
  
  toplot_2 <- as.data.table(p_toplot)[val > floor(num_vecs/2)]
  df_2 <- as.data.table(p_df)
}

ggplot(data = toplot_1, aes(x = x, y = y, color = type)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Point type")
ggsave(paste0("levels_lower_p6m", num_vecs, ".png"))

toplot_1$type <- paste("Old_", toplot_1$type, sep = "")
toplot_2 <- bind_rows(toplot_2, toplot_1)

ggplot(data = toplot_2, aes(x = x, y = y, color = type)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend)) + 
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Point type")
ggsave(paste0("levels_upper_p6m", num_vecs, ".png"))
