# setwd("../../Documents/MSc courses/7850-APC/ProjectNotes/v1/")

library(fastcluster); library(magrittr); library(tidyverse)
# Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"
# timeTaken <- function(pt, sw) {
#   if (is.null(sw[["start_time"]]) & is.null(sw[["end_time"]])) {
#     paste0("Neither start nor end time were collected") %>% return()
#   }else if (is.null(sw[["end_time"]]) & !is.null(sw[["start_time"]])) {
#     paste0("End time was not collected.") %>% return()
#   }else if (!is.null(sw[["end_time"]]) & is.null(sw[["start_time"]])) {
#     paste0("Start time was not collected.") %>% return()
#   }else {
#     z <- difftime(sw[['end_time']], sw[['start_time']], units = "secs") %>% as.double()
#     m <- 60
#     h <- m^2
#     
#     if (z >= h) {
#       hrs <- trunc(z/h)
#       mins <- trunc(z/m - hrs*m)
#       paste0("\nThe ", pt, " process took ", hrs, " hour(s), ", mins, " minute(s), and ", 
#              round(z - hrs*h - mins*m), " second(s).") %>% return()
#     }else if (z < h & z >= m) {
#       mins <- trunc(z/m)
#       paste0("\nThe ", pt, " process took ", mins, " minute(s) and ", round(z - mins*m), " second(s).") %>% return()
#     }else {
#       paste0("\nThe ", pt, " process took ", round(z), " second(s).") %>% return()
#     }  
#   }
# }

file_index <- 9
num_vecs <- 750

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

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


df <- stats::dist(p_a1)
hc <- fastcluster::hclust(df, members = "ave")
plot(hc)
plot(hc, hang = -1)

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
time_taken <- difftime(stopwatch[['end_time']], stopwatch[['start_time']], units = "secs") %>% as.double()




hc <- hclust.vector(USArrests, "cen")
# squared Euclidean distances
hc$height <- hc$height^2
memb <- cutree(hc, k = 1)
