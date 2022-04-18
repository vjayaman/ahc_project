# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/v1/")

library(magrittr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)

num_vecs <- 5 #15
file_index <- 1 #2

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

traceSequences <- function(s_x3, x, r, seq_list) {
  if (nrow(x) > 0) {
    level_of_earlier_merge <- s_x3[variable == "center" & lvl < x$lvl & seqs == x$seqs]$level
    y <- s_x3[level == level_of_earlier_merge]
    mp <- s_x3[level == y$level & variable != "center"]
    seq_list <- append(seq_list, mp[lvl == 1]$value)
    x <- mp[lvl != 1]
    traceSequences(s_x3, x, r, seq_list)
  }else {
    return(tibble(level = r, elements = list(unlist(seq_list))))
  }
}

# from table of 
#   level | x1 | x2 | center | x1_x | x1_y | x2_x | x2_y | center_x | center_y 
# and a table | x | y | of original sequences
# build 
#   level | variable | value | value_x | value_y | seqs | lvl
# 
buildLvlTbl <- function(s_x2, s_a1) {
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