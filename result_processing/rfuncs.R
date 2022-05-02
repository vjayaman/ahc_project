# setwd("C:/Users/vasen/Documents/MSc courses/7850-APC/ProjectNotes/ahc_project")

library(magrittr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)
library(fastcluster)

initialPairs <- function(seq_results) {
  s_a1 <- lapply(seq_results, function(x) {
    strsplit(x, "\t") %>% unlist() %>% as.double()
  }) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
    as_tibble() %>% set_colnames(c("x", "y"))
  return(s_a1)
}

levelCharNumTbl <- function(s_lvlinds, s_outputs) {
  s_clustering <- s_outputs[s_lvlinds+1]
  s_levels <- lapply(s_outputs[s_lvlinds], function(x_i) {
    strsplit(x_i, split = "level ")[[1]][2] %>% as.double()
  }) %>% unlist()
  
  s_x2 <- lapply(s_clustering, function(x_i) {
    x_i %>% gsub("From ", "", .) %>% gsub(" and ", "\t", .) %>% 
      gsub(" to ", "\t", .) %>% 
      strsplit(., "\t") %>% unlist() %>% t() %>% as.data.frame() %>% 
      select(-1)
  }) %>% bind_rows() %>% set_colnames(c("x1", "x2", "center")) %>% 
    add_column(level = s_levels, .before = 1) %>% 
    sepCols(2, .) %>% sepCols(3, .) %>% sepCols(4, .) %>% as_tibble()
  return(s_x2)
}

# orderedPairFormatting <- function(x1, index) {
#   lapply(1:length(x1), function(i) {
#     x2 <- x1[[i]][nchar(x1[[i]]) > 0]
#     x2[index] %>% gsub("\\]|\\[", "", .) %>% strsplit(., ",") %>% unlist() %>% as.double()
#   }) %>% as.data.frame() %>% t() %>% as.data.frame() %>% as_tibble() %>% 
#     rownames_to_column("val") %>% set_colnames(c("val", "x", "y"))
# }

orderedPairs <- function(x1, index) {
  lapply(1:nrow(x1), function(i) {
    df <- x1[i,..index] %>% gsub("\\]|\\[", "", .) %>% strsplit(., ",") %>% unlist() %>% as.double() %>% t() %>% 
      as.data.table() %>% 
      data.table(x1[i,c("proc","level")], .)
  }) %>% bind_rows() %>% set_colnames(c("proc", "level", "x", "y"))
}

sepCols <- function(col_index, df) {
  cx <- colnames(df)[col_index]
  orig_col <- pull(df, cx)
  a1 <- lapply(orig_col, function(x_i) {
    gsub("\\[", "", x_i) %>% gsub("\\]", "", .) %>% strsplit(., ", ") %>% unlist() %>% 
      as.double() %>% t() %>% as.data.frame()
  }) %>% bind_rows() %>% add_column(orig_col, .before = 1) %>% 
    set_colnames(c(cx, paste(cx, c("_x","_y"), sep = "")))
  
  inner_join(df, a1, by = cx) %>% return()  
}


# from table of 
#   level | x1 | x2 | center | x1_x | x1_y | x2_x | x2_y | center_x | center_y 
# and a table | x | y | of original sequences
# build 
#   level | variable | value | value_x | value_y | seqs | lvl
# 
buildLvlTbl <- function(s_x2, s_a1, is_par) {
  sequences0 <- s_a1 %>% 
    add_column(seqs = paste("seq", str_pad(1:nrow(s_a1), 2, pad = "0"), sep = ""))
  sequences0$chr <- paste("[", formatC(sequences0$x, digits = 2, format = "f"), ", ", 
                          formatC(sequences0$y, digits = 2, format = "f"), "]", sep = "")
  sequences <- sequences0 %>% select(-x, -y)
  
  s_x3 <- s_x2 %>% select(x1, x2, center, level) %>% 
    as.data.table() %>% melt.data.table(id.vars = "level") %>% sepCols(3, .) %>% 
    mutate(across(level, as.double)) %>% unique() %>% 
    left_join(., sequences, by = c("value"="chr")) %>% add_column(lvl = 0)
  s_x3[!is.na(seqs)]$lvl <- 1
  s_cnames <- colnames(s_x3)
  
  if (is_par) {
    s_x3 <- s_x2 %>% select(proc, level, x1, x2, center) %>% 
      melt.data.table(id.vars = c("proc", "level")) %>% 
      merge.data.table(., s_x3, by = c("variable", "value", "level")) %>% 
      select(c(all_of(s_cnames), "proc"))
  }
  
  
  todo <- s_x3[is.na(seqs)]$level %>% unique() %>% sort(decreasing = TRUE)
  j <- s_x3$level %>% max()
  
  for (i in 1:length(todo)) {
    isna_inds <- which(is.na(s_x3[level == todo[i]]$seqs)) # print(s_x3[level == todo[i]])
    # if (length(isna_inds) == 3) {print(paste0(i, ": all 3"))}
    # if (length(isna_inds) == 2) {print("two NA")}
    
    if (length(isna_inds) > 1) {
      for (k in 1:length(isna_inds)) {
        i1 <- s_x3[level == todo[i]][isna_inds[k]]
        if (i1$value %in% sequences$chr) {
          s_x3[level == todo[i]][isna_inds[k]]$seqs <- sequences$seqs[sequences$chr==i1$value]
          s_x3[level == todo[i]][isna_inds[k]]$lvl <- max(s_x3[value == i1$value]$lvl)+1
        }else {
          j <- j+1
          s_x3[level == todo[i]][isna_inds[k]]$seqs <- paste("seq", str_pad(j, 2, pad = "0"), sep = "")
          s_x3[level == todo[i]][isna_inds[k]]$lvl <- max(s_x3$lvl)+1
          sequences <- sequences %>% 
            add_row(tibble(seqs = paste("seq", str_pad(j, 2, pad = "0"), sep = ""), 
                           chr = s_x3[level == todo[i]][isna_inds[k]]$value))
        }  
      }
    }else {# print("one NA") # must be the center
      j <- j+1
      s_x3[level == todo[i]][isna_inds]$seqs <- paste("seq", str_pad(j, 2, pad = "0"), sep = "") 
      s_x3[level == todo[i]][isna_inds]$lvl <- max(s_x3$lvl)+1
      sequences <- sequences %>% 
        add_row(tibble(seqs = paste("seq", str_pad(j, 2, pad = "0"), sep = ""), 
                       chr = s_x3[level == todo[i]][isna_inds]$value))
    }
  }
  return(list(s_x3 = s_x3, seq_tbl = sequences))  
}

# generate table of level and absorbed elements into a group
# designed in a sequential way
# will use this for comparing "accuracy" for sequential and parallel groups later
containedElementsTable <- function(s_x3) {
  for (r in sort(unique(s_x3$level), decreasing = TRUE)) {
    # print(r)
    if (r == max(unique(s_x3$level))) {
      accuracy_tbl <- data.table(
        level = r, elements = list(s_x3[level == r & variable != "center"]$value))
    }else {
      merged_pair <- s_x3[level == r & variable != "center"]
      leaves <- merged_pair[lvl == 1]
      if (nrow(leaves) > 0) {
        seq_list <- list(leaves$value)
      }else {
        seq_list <- list()
      }
      # a pair merged at this level that have children
      x1 <- merged_pair[lvl != 1]
      if (nrow(x1) > 0) {
        for (i in 1:nrow(x1)) {
          seq_list <- append(seq_list, traceSequences(s_x3, x1[i,], r, seq_list))
          df <- tibble(level = r, elements = list(unique(unlist(seq_list))))
          accuracy_tbl <- bind_rows(accuracy_tbl, df)
        }  
      }else { # just merging two leaves
        accuracy_tbl <- bind_rows(
          accuracy_tbl, 
          data.table(level = r, elements = list(s_x3[level == r & variable != "center"]$value)))
      }
    }
  }
  return(accuracy_tbl)
}



traceSequences <- function(s_x3, x, r, seq_list) {
  if (nrow(x) > 0) {
    for (j in 1:nrow(x)) {
      z <- x[j,]
      # when was this value formed in a merge?
      # i.e. when it was a center, at a level before the cuurrent one, 
      # and when the sequence ID was the same
      # level_of_earlier_merge <- s_x3[variable == "center" & lvl < z$lvl & seqs == z$seqs]$level
      level_of_earlier_merge <- s_x3[seqs == z$seqs & variable == "center"]$level
      # collect the triple of x1, x2, pair at this earlier level
      y <- s_x3[level == level_of_earlier_merge]
      # what was the merged pair?
      mp <- s_x3[level == unique(y$level) & variable != "center"]
      # add any children from that pair to the list if they were leaves
      seq_list <- append(seq_list, mp[lvl == 1]$value)
      # if they were not leaves, then we want to trace back to find when they 
      # were merged into an internal node -- will keep doint this until we have 
      # added all leaves to the list
      z <- mp[lvl != 1]
      seq_list <- append(seq_list, traceSequences(s_x3, z, r, seq_list))
    }
  }else {
    # print(list(unlist(seq_list)))
    return(list(unlist(seq_list)))
  }
  return(seq_list)
}


# generate table of level and absorbed elements into a group
# designed in a sequential way
# will use this for comparing "accuracy" for sequential and parallel groups later
parallelContainedElementsTable <- function(s_x3, p) {
  dif_levels <- sort(unique(s_x3$level), decreasing = TRUE)
  for (r in dif_levels) {
    # print(r)
    if (r == max(unique(s_x3$level))) {
      accuracy_tbl <- data.table(
        level = r, elements = list(s_x3[level == r & variable != "center"]$value))
    }else {
      merged_pair <- s_x3[level == r & variable != "center"]
      leaves <- merged_pair[lvl == 1]
      if (nrow(leaves) > 0) {
        seq_list <- list(leaves$value)
      }else {
        seq_list <- list()
      }
      # a pair merged at this level that have children
      x1 <- merged_pair[lvl != 1]
      if (nrow(x1) > 0) {
        for (i in 1:nrow(x1)) {
          seq_list <- append(seq_list, parTraceSequences(s_x3, x1[i,], r, seq_list, p))
          df <- tibble(level = r, elements = list(unique(unlist(seq_list))))
          accuracy_tbl <- bind_rows(accuracy_tbl, df)
        }  
      }else { # just merging two leaves
        accuracy_tbl <- bind_rows(
          accuracy_tbl, 
          data.table(level = r, elements = list(s_x3[level == r & variable != "center"]$value)))
      }
    }
  }
  return(accuracy_tbl)
}


parTraceSequences <- function(s_x3, x, r, seq_list, p) {
  if (nrow(x) > 0) {
    for (j in 1:nrow(x)) {
      z <- x[j,]
      # when was this value formed in a merge?
      # i.e. when it was a center, at a level before the cuurrent one, 
      # and when the sequence ID was the same
      # level_of_earlier_merge <- s_x3[variable == "center" & lvl < z$lvl & seqs == z$seqs]$level
      
      # arbitrarily pick the first of these (since we consider it a leaf if proc != given or 0)
      earlier_merge <- s_x3[seqs == z$seqs & variable == "center"][1]
      
      # collect the triple of x1, x2, pair at this earlier level
      y <- s_x3[level == earlier_merge$level & proc == earlier_merge$proc]
      
      # what was the merged pair?
      mp <- y[variable != "center"]
      
      # add any children from that pair to the list if they were leaves or from a different proc
      seq_list <- append(seq_list, mp[lvl == 1 & proc %in% c(0, p)]$value)
      seq_list <- append(seq_list, mp[lvl == 1 & !(proc %in% c(0,p))]$value)
      # if they were not leaves, then we want to trace back to find when they 
      # were merged into an internal node -- will keep doint this until we have 
      # added all leaves to the list
      z <- mp[lvl != 1 & proc %in% c(0, p)]
      seq_list <- append(seq_list, parTraceSequences(s_x3, z, r, seq_list, p))
    }
  }else {
    # print(list(unlist(seq_list)))
    return(list(unlist(seq_list)))
  }
  return(seq_list)
}


orderedPairToDF <- function(df, col_n) {
  pull(df, col_n) %>% gsub("\\]|\\[", "", .) %>% strsplit(., ",") %>% as.data.frame() %>% 
    t() %>% as.data.table() %>% set_colnames(c("x","y")) %>% 
    mutate(across(c(x,y), as.double)) %>% 
    add_column(df, .before = 1)
}





# PLOTTING --------------------------------------------------------------------------
seqPlotsGeneral <- function(num_vecs, s_outputs, plot_level) {
  file_index <- str_pad(as.double(num_vecs), 4, pad = "0")
  
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
  
  # assertthat::assert_that(nrow(s_toplot[proc != 0]) == 0)
  s_d1 <- s_toplot %>% arrange(proc, -level) %>% filter(pt < 2)
  s_d2a <- s_d1[s_d1$pt == 0] %>% rownames_to_column("row_id") %>% select(-pt)
  
  s_d2b <- s_d1[s_d1$pt == 1] %>% 
    set_colnames(c("proc", "level", "xend","yend","pt")) %>% 
    rownames_to_column("row_id")
  
  s_d2 <- inner_join(s_d2a, s_d2b, by = c("proc","level","row_id")) %>% select(-row_id)
  s_df <- bind_rows(s_d2, s_d2) %>% bind_rows(s_d2)
  s_toplot$pt[s_toplot$pt == "1"] <- "Merged pair"
  s_toplot$pt[s_toplot$pt == "2"] <- "Center"

  threshold <- median(sort(s_toplot$level))
  
  toplot_1 <- as.data.table(s_toplot)[level >= threshold]
  df_1 <- as.data.table(s_df)[level >= threshold]
  
  toplot_2 <- as.data.table(s_toplot)
  df_2 <- as.data.table(s_df)
  
  toplot_1$level <- as.character(toplot_1$level)
  toplot_2$level <- as.character(toplot_2$level)
  df_1$level <- as.character(df_1$level)
  df_2$level <- as.character(df_2$level)
  
  if (plot_level == 1) {
    # fig1 <- paste0("figs_and_tbls/levels_lower_p1m", num_vecs, ".png")
    ggplot(data = toplot_1, aes(x = x, y = y, color = level)) + geom_point() + 
      xlim(min(toplot_2$x), max(toplot_2$x)) + 
      ylim(min(toplot_2$y), max(toplot_2$y)) + 
      geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
      ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2 (halfway)")) + 
      labs(color = "Level (root is 0)") + 
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1, byrow = TRUE))
    # ggsave(fig1)  
  }else {
    # fig2 <- paste0("figs_and_tbls/levels_upper_p1m", num_vecs, ".png")
    ggplot(data = toplot_2, aes(x = x, y = y, color = level)) + geom_point() + 
      xlim(min(toplot_2$x), max(toplot_2$x)) + 
      ylim(min(toplot_2$y), max(toplot_2$y)) + 
      geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend, color = level)) + 
      ggtitle(paste0("Naive sequential implementation for ", num_vecs, " vectors and N = 2s")) + 
      labs(color = "Level (root is 0)") + 
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1, byrow = TRUE))
    # ggsave(fig2) 
  }
}

parPlotsGeneral <- function(num_vecs, p_outputs, plot_level) {
  file_index <- str_pad(as.double(num_vecs), 4, pad = "0")
  
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
  
  if (plot_level == 1) {
    # fig1 <- paste0("figs_and_tbls/levels_lower_p6m", num_vecs, ".png")
    ggplot(data = toplot_1, aes(x = x, y = y, color = level)) + geom_point() + 
      xlim(min(toplot_2$x), max(toplot_2$x)) + 
      ylim(min(toplot_2$y), max(toplot_2$y)) + 
      geom_segment(aes(x = df_1$x, y = df_1$y, xend = df_1$xend, yend = df_1$yend)) + 
      ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
      labs(color = "Level (root is 0)") + 
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1, byrow = TRUE))
    # ggsave(fig1)
  }else {
    # fig2 <- paste0("figs_and_tbls/levels_upper_p6m", num_vecs, ".png")
    ggplot(data = toplot_2, aes(x = x, y = y, color = level)) + geom_point() + 
      xlim(min(toplot_2$x), max(toplot_2$x)) + 
      ylim(min(toplot_2$y), max(toplot_2$y)) + 
      geom_segment(aes(x = df_2$x, y = df_2$y, xend = df_2$xend, yend = df_2$yend)) + 
      ggtitle(paste0("Parallel implementation 1 for ", num_vecs, " vectors and N = 2")) + 
      labs(color = "Level (root is 0)") + 
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1, byrow = TRUE))
    # ggsave(fig2)
  }
}





# TIMING --------------------------------------------------------------------------
parallelFixNVaryP <- function(pfiles, m_val) {
  pfiles <- grep(m_val, pfiles, value = TRUE)
  parallel_times <- lapply(1:length(pfiles), function(i) {
    x <- pfiles[i]
    p_outputs <- readLines(x)
    p_ind1 <- which(grepl("Total time", p_outputs))
    if (length(p_ind1) > 0) {
      y <- strsplit(names(x), split = "m") %>% unlist() %>% gsub("p", "", .) %>% as.double()
      z <- strsplit(p_outputs[p_ind1], " ")[[1]][4] %>% as.double()
      data.frame(t(y), z) %>% set_colnames(c("num_nodes", "num_vecs", "time")) 
    }
  }) %>% bind_rows()
  return(parallel_times)
}

runParallelFixNVaryP <- function(fpath) {
  pfiles <- list.files(fpath, full.names = TRUE)
  pnames <- pfiles %>% gsub(fpath, "", .) %>% gsub(".txt", "", .)
  names(pfiles) <- pnames
  
  parallelFixNVaryP(pfiles, "m015") %>% 
    bind_rows(., parallelFixNVaryP(pfiles, "m060")) %>% 
    bind_rows(., parallelFixNVaryP(pfiles, "m300")) %>% as.data.table() %>% return()
}

seqFixNVaryP <- function(sfiles, m_val) {
  sfiles <- grep(m_val, sfiles, value = TRUE)
  lapply(1:length(sfiles), function(i) {
    x <- sfiles[i]
    s_outputs <- readLines(x)
    s_ind1 <- which(grepl("Total time", s_outputs))
    y <- strsplit(names(x), split = "m") %>% unlist() %>% gsub("p", "", .) %>% as.double()
    z <- strsplit(s_outputs[s_ind1], " ")[[1]][4] %>% as.double()
    data.frame(num_nodes = y[1], num_vecs = y[2], time = z)
  }) %>% bind_rows() %>% return()
  # return(naive_sequential_times)
}

runseqFixNVaryP <- function(fpath) {
  sfiles <- list.files(fpath, full.names = TRUE, pattern = ".txt")
  snames <- sfiles %>% gsub(fpath, "", .) %>% gsub(".txt", "", .)
  names(sfiles) <- snames
  # naive_sequential_times <- 
  seqFixNVaryP(sfiles, "m015") %>% 
    bind_rows(., seqFixNVaryP(sfiles, "m060")) %>% 
    bind_rows(., seqFixNVaryP(sfiles, "m300")) %>% as.data.table() %>% return()
}

runBestSeq <- function(fpath) {
  num_vec_list <- c(15, 60, 300)
  # sequential_times <- 
  lapply(1:length(num_vec_list), function(i) {
    file_index <- i
    num_vecs <- num_vec_list[i]
    
    start_time <- Sys.time()
    
    sfiles <- list.files(fpath, full.names = TRUE, pattern = ".txt")
    snames <- sfiles %>% gsub(fpath, "", .) %>% gsub(".txt", "", .)
    names(sfiles) <- snames
    s_outputs <- readLines(sfiles[grep(as.character(num_vecs), snames)])
    
    s_ind1 <- which(grepl("Matrix A: ", s_outputs))+1
    seq_results <- s_outputs[s_ind1:(s_ind1+(num_vecs-1))]
    
    s_a1 <- lapply(seq_results, function(x) {
      strsplit(x, "\t") %>% unlist() %>% as.double()
    }) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
      as_tibble() %>% set_colnames(c("x", "y"))
    
    df <- stats::dist(s_a1)
    hc <- fastcluster::hclust(df, members = "ave") # plot(hc); plot(hc, hang = -1)
    
    end_time <- Sys.time() - start_time
    
    data.frame(num_nodes = 1, num_vecs, time = as.character(end_time) %>% as.double())
  }) %>% bind_rows() %>% as.data.table() %>% return()
}





