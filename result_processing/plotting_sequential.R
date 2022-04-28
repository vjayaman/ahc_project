# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/ahc_project/")

source("result_processing/rfuncs.R")

num_vecs <- 15
file_index <- 3

# ------------------------------------------------------------------------------
# SEQUENTIAL -------------------------------------------------------------------

sfiles <- list.files("outputs/seq/", full.names = TRUE)
snames <- sfiles %>% gsub("outputs/seq/", "", .) %>% gsub(".txt", "", .)
names(sfiles) <- snames
s_outputs <- readLines(sfiles[file_index])

s_ind1 <- which(grepl("Matrix A: ", s_outputs))+1
seq_results <- s_outputs[s_ind1:(s_ind1+(num_vecs-1))]

s_a1 <- initialPairs(seq_results)

s_pairs <- grep("From ", s_outputs, value = TRUE)
s_clustering <- grep("Minimum distance ", s_outputs, value = TRUE)

# note that the "levels" are decreasing --> max is at the leaves, 1 is at the root

s_levels <- lapply(s_clustering, function(x_i) {
  t_i <- strsplit(x_i, split = "level ")[[1]][2]
  u_i <- strsplit(t_i, ", process ")[[1]][1] %>% as.double()
  v_i <- strsplit(strsplit(t_i, ", process ")[[1]][2], "\t")[[1]][1] %>% as.double()
  data.table(level = u_i, proc = v_i)
}) %>% bind_rows()


s_x2 <- lapply(s_pairs, function(x_i) {
  y_i <- strsplit(x_i, "From ")[[1]][2] %>% 
    # gsub("From ", "", .) %>% 
    gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
  y_i[nchar(y_i) > 0] %>% t() %>% as.data.frame()
}) %>% bind_rows() %>% 
  set_colnames(c("x1", "x2", "center")) %>% 
  add_column(level = s_levels$level, .before = 1) %>% 
  add_column(proc = s_levels$proc, .before = 1) %>% 
  sepCols(3, .) %>% sepCols(4, .) %>% sepCols(5, .) %>% as_tibble() %>% unique() %>% 
  as.data.table()

s_x2[proc != 0]$level <- s_x2[proc != 0]$level + max(s_x2$level) - 1
s_x2$level <- s_x2$level - 2 # bring to range (0,max)

s_first_pair <- s_x2 %>% select(proc, level, x1) %>% orderedPairToDF(., "x1") %>% 
  select(proc, level, x, y) %>% add_column(type = 1)
s_second_pair <- s_x2 %>% select(proc, level, x2) %>% orderedPairToDF(., "x2") %>% 
  select(proc, level, x, y) %>% add_column(type = 1)
s_new_center <- s_x2 %>% select(proc, level, center) %>% orderedPairToDF(., "center") %>% 
  select(proc, level, x, y) %>% add_column(type = 2)


s_toplot <- bind_rows(s_first_pair, s_second_pair) %>% bind_rows(., s_new_center) %>% 
  as.data.table() %>% mutate(across(level, as.double))


s_d1 <- s_toplot[type == 1] %>% arrange(level, proc) %>% 
  add_column(new_val = rep(c(0,1), nrow(s_toplot[s_toplot$type == 1,])/2))

s_d2a <- s_d1[new_val == 0, c("proc","level","x","y","type")]
s_d2b <- s_d1[new_val == 1, c("proc","level","x","y","type")] %>% set_colnames(c("proc","level","xend","yend","type"))
s_df <- inner_join(s_d2a, s_d2b, by = c("proc","level","type")) %>% 
  mutate(across(level, as.double)) # %>% unique()
s_df <- bind_rows(s_df, s_df) %>% bind_rows(., s_df) # geom_segment requires same number of rows

s_toplot$type[s_toplot$type == "1"] <- "Merged pair"
s_toplot$type[s_toplot$type == "2"] <- "Center"

threshold <- median(unique(s_toplot$level))

toplot_1 <- as.data.table(s_toplot) %>% filter(level > threshold)
df_1 <- as.data.table(s_df)[level > threshold]

toplot_2 <- as.data.table(s_toplot)[level <= threshold]
df_2 <- as.data.table(s_df)

toplot_1$level <- as.character(toplot_1$level)
toplot_2$level <- as.character(toplot_2$level)
toplot_2 <- bind_rows(toplot_2, toplot_1)

ggplot(data = toplot_1, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
  ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2 (halfway)")) + 
  labs(color = "Level (root is 0)")
ggsave(paste0("figs_and_tbls/levels_lower_p1m", num_vecs, ".png"))



ggplot(data = toplot_2, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend, color = level)) + 
  ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2s")) + 
  labs(color = "Level (root is 0)")
ggsave(paste0("figs_and_tbls/levels_upper_p1m", num_vecs, ".png"))











