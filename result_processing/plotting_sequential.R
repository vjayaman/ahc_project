# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/ahc_project/")

source("result_processing/rfuncs.R")

num_vecs <- 15
file_index <- str_pad(num_vecs, 4, pad = "0")

sfiles <- list.files("outputs/seq/", full.names = TRUE)
snames <- sfiles %>% gsub("outputs/seq/", "", .) %>% gsub(".txt", "", .)
names(sfiles) <- snames
s_outputs <- readLines(sfiles[grep(file_index, snames)])

s_ind1 <- which(grepl("Matrix A: ", s_outputs))+1
seq_results <- s_outputs[s_ind1:(s_ind1+(num_vecs-1))]

s_a1 <- initialPairs(seq_results)

s_pairs <- grep("From ", s_outputs, value = TRUE)
s_clustering <- grep("Minimum distance ", s_outputs, value = TRUE)

# note that the "levels" are decreasing --> max is at the leaves, 1 is at the root

s_levels <- lapply(s_clustering, function(x_i) {
  y_i <- strsplit(strsplit(x_i, split = "level ")[[1]][2], ", process ")[[1]]
  y_i_proc <- strsplit(y_i[2], split = ": ")[[1]][1] %>% as.double()
  data.table(level = as.double(y_i[1]), proc = y_i_proc)
}) %>% bind_rows()

df <- lapply(s_pairs, function(x_i) {
  y_i <- strsplit(x_i, split = ": ")[[1]][2] %>% 
    gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
  y_i[nchar(y_i) > 0] %>% t() %>% as.data.frame()
}) %>% bind_rows()

s_x2 <- df %>% set_colnames(c("x1", "x2", "center")) %>% 
  add_column(level = s_levels$level, .before = 1) %>% 
  add_column(proc = s_levels$proc, .before = 1) %>% 
  sepCols(3, .) %>% sepCols(4, .) %>% sepCols(5, .) %>% as_tibble() %>% unique() %>% 
  as.data.table()

s_first_pair <- orderedPairs(s_x2, "x1") %>% add_column(pt = "0")
s_second_pair <- orderedPairs(s_x2, "x2") %>% add_column(pt = "1")
s_new_center <- orderedPairs(s_x2, "center") %>% add_column(pt = "2")

s_toplot <- bind_rows(s_first_pair, s_second_pair) %>% 
  bind_rows(., s_new_center) %>% 
  mutate(level = level - min(level))

assertthat::assert_that(nrow(s_toplot[proc != 0]) == 0)

s_d1 <- s_toplot %>% arrange(proc, -level) %>% filter(pt < 2)


s_d2a <- s_d1[s_d1$pt == 0] %>% rownames_to_column("row_id") %>% select(-pt)

s_d2b <- s_d1[s_d1$pt == 1] %>% 
  set_colnames(c("proc", "level", "xend","yend","pt")) %>% 
  rownames_to_column("row_id")

s_d2 <- inner_join(s_d2a, s_d2b, by = c("proc","level","row_id")) %>% select(-row_id)

s_df <- bind_rows(s_d2, s_d2) %>% bind_rows(s_d2)

s_toplot$pt[s_toplot$pt == "1"] <- "Merged pair"
s_toplot$pt[s_toplot$pt == "2"] <- "Center"

s_toplot$type[s_toplot$type == "1"] <- "Merged pair"
s_toplot$type[s_toplot$type == "2"] <- "Center"


threshold <- median(sort(s_toplot$level))

toplot_1 <- as.data.table(s_toplot)[level >= threshold]
df_1 <- as.data.table(s_df)[level >= threshold]

toplot_2 <- as.data.table(s_toplot)
df_2 <- as.data.table(s_df)

toplot_1$level <- as.character(toplot_1$level)
toplot_2$level <- as.character(toplot_2$level)
df_1$level <- as.character(df_1$level)
df_2$level <- as.character(df_2$level)


ggplot(data = toplot_1, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
  ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2 (halfway)")) + 
  labs(color = "Level (root is 0)")
# ggsave(paste0("figs_and_tbls/levels_lower_p1m", num_vecs, ".png"))



ggplot(data = toplot_2, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend, color = level)) + 
  ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2s")) + 
  labs(color = "Level (root is 0)")
# ggsave(paste0("figs_and_tbls/levels_upper_p1m", num_vecs, ".png"))











