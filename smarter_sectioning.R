s_a1$color <- "a"
mx <- mean(s_a1$x)
my <- mean(s_a1$y)
for (i in 1:nrow(s_a1)) {
  if (s_a1$x[i] <= mx & s_a1$y[i] <= my) {
    s_a1$color[i] <- "b"
  }else if (s_a1$x[i] <= mx & s_a1$y[i] > my) {
    s_a1$color[i] <- "c"
  }else if (s_a1$x[i] > mx & s_a1$y[i] > my) {
    s_a1$color[i] <- "d"
  }
}

ggplot(s_a1, aes(x, y, color = color)) + geom_point()

