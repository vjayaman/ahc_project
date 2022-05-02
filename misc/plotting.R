# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/v1/")

library(magrittr)
library(tidyverse)
library(ggplot2)

num_vecs <- 15
file_index <- 2

orderedPairFormatting <- function(x1, index) {
  lapply(1:length(x1), function(i) {
    x1[[i]][index] %>% gsub("\\]|\\[", "", .) %>% strsplit(., ",") %>% unlist() %>% as.double()
  }) %>% as.data.frame() %>% t() %>% as.data.frame() %>% as_tibble() %>% 
    rownames_to_column("val") %>% set_colnames(c("val", "x", "y"))
}
# ---------------Third attempt ----------------------------------


# PARALLEL
pfiles <- list.files("outputs/parallel/", full.names = TRUE)
pnames <- pfiles %>% gsub("outputs/parallel/", "", .) %>% gsub(".txt", "", .)
names(pfiles) <- pnames
p_outputs <- readLines(pfiles[file_index])

# for (x in pfiles) {
#   p_outputs <- readLines(x)
#   p_ind1 <- which(grepl("Total time", p_outputs))
#   print(p_outputs[p_ind1])
# }

p_ind1 <- which(grepl("Matrix A: ", p_outputs))[1]+1
p_matA <- p_outputs[p_ind1:(p_ind1+(num_vecs-1))]
p_a1 <- lapply(p_matA, function(x) {
    strsplit(x, "\t") %>% unlist() %>% as.double()
  }) %>% as.data.frame() %>% t() %>% as.data.frame() %>% 
    as_tibble() %>% set_colnames(c("x", "y"))

p_clustering <- p_outputs[which(grepl("From ", p_outputs))]
p_x1 <- lapply(p_clustering, function(x_i) {
  x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
})


p_first_pair <- orderedPairFormatting(p_x1, 2) %>% add_column(type = "1")
p_second_pair <- orderedPairFormatting(p_x1, 3) %>% add_column(type = "1")
p_new_center <- orderedPairFormatting(p_x1, 4) %>% add_column(type = "2")

p_toplot <- bind_rows(p_first_pair, p_second_pair) %>% bind_rows(., p_new_center)

p_d1 <- p_toplot[p_toplot$type == 1,] %>% arrange(val) %>% 
  add_column(new_val = rep(c(0,1), nrow(p_toplot[p_toplot$type == 1,])/2))

p_d2a <- p_d1[p_d1$new_val == 0,c("val","x","y","type")]
p_d2b <- p_d1[p_d1$new_val == 1,c("val","x","y","type")] %>% 
  set_colnames(c("val","xend","yend","type"))
p_d2 <- inner_join(p_d2a, p_d2b, by = c("val","type"))

p_df <- bind_rows(p_d2, p_d2) %>% bind_rows(p_d2)

p_toplot$type[p_toplot$type == "1"] <- "Merged pair"
p_toplot$type[p_toplot$type == "2"] <- "Center"
p_plot <- 
  ggplot(data = p_toplot, aes(x = x, y = y, color = type)) + geom_point() + 
  geom_segment(aes(x = p_df$x, y = p_df$y, xend = p_df$xend, yend = p_df$yend)) + 
  ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Point type")
p_plot
ggsave(paste0("parallel_plot.png"))


# i <- 72
#   ggplot(data = p_toplot, aes(x = x, y = y, color = type)) + geom_point() + 
#   geom_segment(aes(x = p_df$x[i], y = p_df$y[i], xend = p_df$xend[i], yend = p_df$yend[i])) + 
#   ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
#   labs(color = "Point type")


# SEQUENTIAL

sfiles <- list.files("outputs/seq/", full.names = TRUE)
snames <- sfiles %>% gsub("outputs/seq/", "", .) %>% gsub(".txt", "", .)
names(sfiles) <- snames
s_outputs <- readLines(sfiles[file_index])

s_ind1 <- which(grepl("Matrix A: ", s_outputs))+1
seq_results <- s_outputs[s_ind1:(s_ind1+(num_vecs-1))]

s_a1 <- lapply(seq_results, function(x) {
  strsplit(x, "\t") %>% unlist() %>% as.double()
}) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
  as_tibble() %>% set_colnames(c("x", "y"))

s_clustering <- s_outputs[which(grepl("From ", s_outputs))]
s_x1 <- lapply(s_clustering, function(x_i) {
  x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
    gsub(" to ", "\t", .) %>% 
    strsplit(., "\t") %>% unlist()
})


s_first_pair <- orderedPairFormatting(s_x1, 2) %>% add_column(type = "1")
s_second_pair <- orderedPairFormatting(s_x1, 3) %>% add_column(type = "1")
s_new_center <- orderedPairFormatting(s_x1, 4) %>% add_column(type = "2")

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
s_plot <- ggplot(data = s_toplot, aes(x = x, y = y, color = type)) + geom_point() + 
  geom_segment(aes(x = s_df$x, y = s_df$y, xend = s_df$xend, yend = s_df$yend)) + 
  ggtitle(paste0("Sequential implementation 1 for ", num_vecs, " vectors and N = 2")) + 
  labs(color = "Point type")
s_plot
ggsave(paste0("sequential_plot.png"))


# for (i in 1:length(unique(toplot$val))) {
#   df <- d2[d2$val == i,]
#   # df_line <- toplot[toplot$val <= i & toplot$type == 1,]
#   p <- p +
#       geom_point(data = toplot[toplot$val == i,], aes (x = x, y = y, color = type)) + 
#       geom_segment(aes(x = df$x, y = df$y, xend = df$xend, yend = df$yend))
#     p
#     # y = toplot[toplot$val == i & toplot$type == 1,]$y[1], 
#     # ggplot(a1, aes(x = x, y = y)) + geom_point() + 
#     #   geom_point(data = toplot[toplot$val < i,], aes (x = x, y = y, color = "red")) + 
#       
#   ggsave(paste0("plot_", i, ".png"))
# }
