# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/ahc_project/")

source("result_processing/rfuncs.R")

num_vecs <- 15
file_index <- str_pad(num_vecs, 4, pad = "0")

pfiles <- list.files("outputs/parallel/", full.names = TRUE)
pnames <- pfiles %>% gsub("outputs/parallel/", "", .) %>% gsub(".txt", "", .)
names(pfiles) <- pnames
p_outputs <- readLines(pfiles[grep(file_index, pnames)])

p_ind1 <- which(grepl("Matrix A: ", p_outputs))[1]+1
p_matA <- p_outputs[p_ind1:(p_ind1+(num_vecs-1))]

p_a1 <- lapply(p_matA, function(x) {
  strsplit(x, "\t") %>% unlist() %>% as.double()
}) %>% as.data.frame() %>% t() %>% as.data.frame() %>% 
  as_tibble() %>% set_colnames(c("x", "y"))


p_pairs <- grep("From ", p_outputs, value = TRUE)
# p_pairs <- gsub("\t", ": ", p_pairs)

p_clustering <- grep("Minimum distance ", p_outputs, value = TRUE)
# p_clustering <- gsub("\t", ": ", p_clustering)

# note that the "levels" are decreasing --> max is at the leaves, 1 is at the root

p_levels <- lapply(p_clustering, function(x_i) {
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
  add_column(level = p_levels$level, .before = 1) %>% 
  add_column(proc = p_levels$proc, .before = 1) %>% 
  sepCols(3, .) %>% sepCols(4, .) %>% sepCols(5, .) %>% as_tibble() %>% unique() %>% 
  as.data.table()

p_first_pair <- orderedPairs(p_x2, "x1") %>% add_column(pt = "0")
p_second_pair <- orderedPairs(p_x2, "x2") %>% add_column(pt = "1")
p_new_center <- orderedPairs(p_x2, "center") %>% add_column(pt = "2")

p_toplot <- bind_rows(p_first_pair, p_second_pair) %>% 
  bind_rows(., p_new_center) %>% 
  mutate(level = level - min(level))

p_toplot[proc != 0] <- p_toplot[proc != 0] %>% 
  mutate(level = level + max(p_toplot[proc==0]$level) + 1)

p_d1 <- p_toplot %>% arrange(proc, -level) %>% filter(pt < 2)


p_d2a <- p_d1[p_d1$pt == 0] %>% rownames_to_column("row_id") %>% select(-pt)

p_d2b <- p_d1[p_d1$pt == 1] %>% 
  set_colnames(c("proc", "level", "xend","yend","pt")) %>% 
  rownames_to_column("row_id")

p_d2 <- inner_join(p_d2a, p_d2b, by = c("proc","level","row_id")) %>% select(-row_id)

p_df <- bind_rows(p_d2, p_d2) %>% bind_rows(p_d2)

p_toplot$pt[p_toplot$pt == "1"] <- "Merged pair"
p_toplot$pt[p_toplot$pt == "2"] <- "Center"

threshold <- median(sort(p_toplot$level))

toplot_1 <- as.data.table(p_toplot)[level >= threshold]
df_1 <- as.data.table(p_df)[level >= threshold]

toplot_2 <- as.data.table(p_toplot)
df_2 <- as.data.table(p_df)

toplot_1$level <- as.character(toplot_1$level)
toplot_2$level <- as.character(toplot_2$level)
df_1$level <- as.character(df_1$level)
df_2$level <- as.character(df_2$level)

ggplot(data = toplot_1, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Level (root is 0)") + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
# ggsave(paste0("levels_lower_p6m", num_vecs, ".png"))


ggplot(data = toplot_2, aes(x = x, y = y, color = level)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend)) + 
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Level (root is 0)") + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
# ggsave(paste0("levels_upper_p6m", num_vecs, ".png"))
