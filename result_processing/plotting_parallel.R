# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/ahc_project/")

source("result_processing/rfuncs.R")

<<<<<<< HEAD
num_vecs <- 5
file_index <- 1
=======
num_vecs <- 15
file_index <- 3
>>>>>>> cd467e90d28f2299ac04a6b9729378976d9cc408


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


p_pairs <- grep("From ", p_outputs, value = TRUE)
p_clustering <- grep("Minimum distance ", p_outputs, value = TRUE)

# note that the "levels" are decreasing --> max is at the leaves, 1 is at the root
p_levels <- lapply(p_clustering, function(x_i) {
<<<<<<< HEAD
  y_i <- strsplit(strsplit(x_i, split = "level ")[[1]][2], ", process ")[[1]]
  y_i_proc <- strsplit(y_i[2], split = ": ")[[1]][1] %>% as.double()
  data.table(level = as.double(y_i[1]), proc = y_i_proc)
}) %>% bind_rows()

df <- lapply(p_pairs, function(x_i) {
  y_i <- strsplit(x_i, split = ": ")[[1]][2] %>% 
    gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
  y_i[nchar(y_i) > 0] %>% t() %>% as.data.frame()
}) %>% bind_rows() 

p_x2 <- df %>% set_colnames(c("x1", "x2", "center")) %>% 
=======
  t_i <- strsplit(x_i, split = "level ")[[1]][2]
  u_i <- strsplit(t_i, ", process ")[[1]][1] %>% as.double()
  v_i <- strsplit(strsplit(t_i, ", process ")[[1]][2], "\t")[[1]][1] %>% as.double()
  data.table(level = u_i, proc = v_i)
}) %>% bind_rows()


p_x2 <- lapply(p_pairs, function(x_i) {
  y_i <- strsplit(x_i, "From ")[[1]][2] %>% 
    # gsub("From ", "", .) %>% 
    gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
  y_i[nchar(y_i) > 0] %>% t() %>% as.data.frame()
}) %>% bind_rows() %>% 
  set_colnames(c("x1", "x2", "center")) %>% 
>>>>>>> cd467e90d28f2299ac04a6b9729378976d9cc408
  add_column(level = p_levels$level, .before = 1) %>% 
  add_column(proc = p_levels$proc, .before = 1) %>% 
  sepCols(3, .) %>% sepCols(4, .) %>% sepCols(5, .) %>% as_tibble() %>% unique() %>% 
  as.data.table()

p_x2[proc != 0]$level <- p_x2[proc != 0]$level + max(p_x2$level) - 1
p_x2$level <- p_x2$level - 2 # bring to range (0,max)

<<<<<<< HEAD
p_x1 <- lapply(p_pairs, function(x_i) {
  strsplit(x_i, split = ": ")[[1]][2] %>% 
    gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
})

# p_first_pair <- orderedPairFormatting(p_x1, 1) %>% add_column(type = "1")
# p_second_pair <- orderedPairFormatting(p_x1, 2) %>% add_column(type = "1")
# p_new_center <- orderedPairFormatting(p_x1, 3) %>% add_column(type = "2")

p_first_pair <- orderedPairs(p_x2, "x1") %>% add_column(pt = "1")
p_second_pair <- orderedPairs(p_x2, "x2") %>% add_column(pt = "1")
p_new_center <- orderedPairs(p_x2, "center") %>% add_column(pt = "2")


p_toplot <- bind_rows(p_first_pair, p_second_pair) %>% bind_rows(., p_new_center) %>% 
  mutate(level = level - 2)

p_d1 <- p_toplot[pt == 1] %>% arrange(proc, -level)
p_d1 <- p_d1 %>% add_column(new_val = rep(c(0,1), nrow(p_d1)/2))

# p_d1 <- p_toplot[p_toplot$type == 1,] %>% arrange(val) %>% 
#   add_column(new_val = rep(c(0,1), nrow(p_toplot[p_toplot$type == 1,])/2))
# p_d2a <- p_d1[p_d1$new_val == 0,c("val","x","y","type")]
# p_d2b <- p_d1[p_d1$new_val == 1,c("val","x","y","type")] %>% 
#   set_colnames(c("val","xend","yend","type"))

p_d2a <- p_d1[p_d1$new_val == 0] %>% select(-new_val)

p_d2b <- p_d1[p_d1$new_val == 1] %>% select(-new_val) %>% 
  set_colnames(c("proc", "level", "xend","yend","pt"))

p_d2 <- inner_join(p_d2a, p_d2b, by = c("proc","level","pt"))

p_df <- bind_rows(p_d2, p_d2) %>% bind_rows(p_d2)
=======
p_first_pair <- p_x2 %>% select(proc, level, x1) %>% orderedPairToDF(., "x1") %>% 
  select(proc, level, x, y) %>% add_column(type = 1)
p_second_pair <- p_x2 %>% select(proc, level, x2) %>% orderedPairToDF(., "x2") %>% 
  select(proc, level, x, y) %>% add_column(type = 1)
p_new_center <- p_x2 %>% select(proc, level, center) %>% orderedPairToDF(., "center") %>% 
  select(proc, level, x, y) %>% add_column(type = 2)
  

p_toplot <- bind_rows(p_first_pair, p_second_pair) %>% bind_rows(., p_new_center) %>% 
  as.data.table() %>% mutate(across(level, as.double))


p_d1 <- p_toplot[type == 1] %>% arrange(level, proc) %>% 
  add_column(new_val = rep(c(0,1), nrow(p_toplot[p_toplot$type == 1,])/2))

p_d2a <- p_d1[new_val == 0, c("proc","level","x","y","type")]
p_d2b <- p_d1[new_val == 1, c("proc","level","x","y","type")] %>% set_colnames(c("proc","level","xend","yend","type"))
p_df <- inner_join(p_d2a, p_d2b, by = c("proc","level","type")) %>% 
  mutate(across(level, as.double)) # %>% unique()
p_df <- bind_rows(p_df, p_df) %>% bind_rows(., p_df) # geom_segment requires same number of rows
>>>>>>> cd467e90d28f2299ac04a6b9729378976d9cc408

p_toplot$type[p_toplot$pt == "1"] <- "Merged pair"
p_toplot$type[p_toplot$pt == "2"] <- "Center"


<<<<<<< HEAD
# p_toplot <- p_toplot %>% mutate(across(val, as.double))
# p_df <- p_df %>% mutate(across(val, as.double))

threshold <- median(sort(unique(p_toplot$level)))

toplot_1 <- as.data.table(p_toplot)[level >= threshold]
df_1 <- as.data.table(p_df)[level >= threshold]

toplot_2 <- as.data.table(p_toplot)
=======
threshold <- median(unique(p_toplot$level))

toplot_1 <- as.data.table(p_toplot) %>% filter(level > threshold)
df_1 <- as.data.table(p_df)[level > threshold]

toplot_2 <- as.data.table(p_toplot)[level <= threshold]
>>>>>>> cd467e90d28f2299ac04a6b9729378976d9cc408
df_2 <- as.data.table(p_df)

toplot_1$level <- as.character(toplot_1$level)
toplot_2$level <- as.character(toplot_2$level)
<<<<<<< HEAD
df_1$level <- as.character(df_1$level)
df_2$level <- as.character(df_2$level)
=======
toplot_2 <- bind_rows(toplot_2, toplot_1)

>>>>>>> cd467e90d28f2299ac04a6b9729378976d9cc408

ggplot(data = toplot_1, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
<<<<<<< HEAD
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Level (root is 0)")
# ggsave(paste0("levels_lower_p6m", num_vecs, ".png"))


=======
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2 (halfway)")) + 
  labs(color = "Level (root is 0)")
ggsave(paste0("figs_and_tbls/levels_lower_p6m", num_vecs, ".png"))



>>>>>>> cd467e90d28f2299ac04a6b9729378976d9cc408
ggplot(data = toplot_2, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend, color = level)) + 
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Level (root is 0)")
<<<<<<< HEAD
# ggsave(paste0("levels_upper_p6m", num_vecs, ".png"))
=======
ggsave(paste0("figs_and_tbls/levels_upper_p6m", num_vecs, ".png"))



# ggplot(p_toplot, aes(x = x, y = y, color = type)) + geom_point() + 
#   geom_segment(aes(x = p_df[val <= 1]$x, y = p_df[val <= 1]$y, 
#                    xend = p_df[val <= 1]$xend, yend = p_df[val <= 1]$yend)) + 
#   geom_segment(aes(x = p_df[val == 2]$x, y = p_df[val == 2]$y, 
#                    xend = p_df[val == 2]$xend, yend = p_df[val == 2]$yend))

# for (x in sort(unique(p_toplot$val))) {
#   a1 <- p_toplot[which(p_toplot$val <= x)]
#   df_v <- p_df %>% select(-type) %>% merge.data.table(., a1, by = c("val", "x", "y"))
#   ggplot(data = df_v, aes(x = x, y = y, color = type)) + geom_point() + 
#     xlim(min(p_toplot$x), max(p_toplot$x)) +
#     ylim(min(p_toplot$y), max(p_toplot$y)) +
#     geom_segment(aes(x = x, y = y,xend = xend, yend = yend)) +
#     ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
#     labs(color = "Point type")
#   ggsave(paste0("plot_", x, "_", num_vecs, ".png"))
# }

>>>>>>> cd467e90d28f2299ac04a6b9729378976d9cc408
