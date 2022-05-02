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

s_clustering <- s_outputs[which(grepl("From ", s_outputs))]
s_x1 <- lapply(s_clustering, function(x_i) {
  x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
})


s_first_pair <- orderedPairFormatting(s_x1, 1) %>% add_column(type = "1")
s_second_pair <- orderedPairFormatting(s_x1, 2) %>% add_column(type = "1")
s_new_center <- orderedPairFormatting(s_x1, 3) %>% add_column(type = "2")

s_toplot <- bind_rows(s_first_pair, s_second_pair) %>% bind_rows(., s_new_center)

s_d1 <- s_toplot[s_toplot$type == 1,] %>% arrange(val) %>% 
  add_column(new_val = rep(c(0,1), nrow(s_toplot[s_toplot$type == 1,])/2))

s_d2a <- s_d1[s_d1$new_val == 0,c("val","x","y","type")]
s_d2b <- s_d1[s_d1$new_val == 1,c("val","x","y","type")] %>% 
  set_colnames(c("val","xend","yend","type"))
s_d2 <- inner_join(s_d2a, s_d2b, by = c("val","type"))

s_df <- bind_rows(s_d2, s_d2) %>% bind_rows(s_d2)

s_toplot$type[s_toplot$type == "1"] <- "Merged pair"
s_toplot$type[s_toplot$type == "2"] <- "Center"

s_toplot <- s_toplot %>% mutate(across(val, as.double))
s_df <- s_df %>% mutate(across(val, as.double))

toplot_1 <- as.data.table(s_toplot)[val <= floor(num_vecs/2)]
df_1 <- as.data.table(s_df)[val <= floor(num_vecs/2)]

toplot_2 <- as.data.table(s_toplot)
df_2 <- as.data.table(s_df)

ggplot(data = toplot_1, aes(x = x, y = y, color = type)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) + 
  ylim(min(toplot_2$y), max(toplot_2$y)) + 
  geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
  ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Point type")
ggsave(paste0("levels_lower_p1m", num_vecs, ".png"))


ggplot(data = toplot_2, aes(x = x, y = y, color = type)) + geom_point() + 
  xlim(min(toplot_2$x), max(toplot_2$x)) +
  ylim(min(toplot_2$y), max(toplot_2$y)) +
  geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend)) + 
  ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Point type")
ggsave(paste0("levels_upper_p1m", num_vecs, ".png"))















