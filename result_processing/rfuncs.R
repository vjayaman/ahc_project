# setwd("C:/Users/vasen/Documents/MSc courses/7850-APC/ProjectNotes/ahc_project")

library(magrittr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)


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

orderedPairFormatting <- function(x1, index) {
  lapply(1:length(x1), function(i) {
    x2 <- x1[[i]][nchar(x1[[i]]) > 0]
    x2[index] %>% gsub("\\]|\\[", "", .) %>% strsplit(., ",") %>% unlist() %>% as.double()
  }) %>% as.data.frame() %>% t() %>% as.data.frame() %>% as_tibble() %>% 
    rownames_to_column("val") %>% set_colnames(c("val", "x", "y"))
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

