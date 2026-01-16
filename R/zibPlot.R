# method for plot() that plots the distribution

plot.zibSim <- function(x, microbe_name = "Microbe", ...) {

  # get the data
  df <- x$rel_abundance

  # create the plot
  # par(new = TRUE)
  breaks <- seq(0, 1, by = 0.025)
  non_zero <- subset(df, df != 0)

  fit <- MASS::fitdistr(non_zero, dbeta, start = list(shape1 = 1, shape2 = 1))

  zero_prop <- mean(df == 0)
  par(mar = c(5, 4, 4, 5) + 0.1)
  par(col.axis = "blue", col.lab = "blue")
  hist(df, freq = FALSE, xlim = c(0, 1), breaks = breaks,
       main = microbe_name, xlab = "Abundance", ylab = "Density")
  curve(dbeta(x, fit$estimate["shape1"], fit$estimate["shape2"]),
        col = "blue", lwd = 2, add = TRUE)
  par(new = TRUE, col.axis = "red")
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = range(df), ylim = c(0, 1))
  lines(c(0, 0), c(0, zero_prop), col = "red", lwd = 3)
  axis(side = 4, at = seq(0, 1, 0.2))
  mtext("Zero proportion", side = 4, line = 3, col = "red")
  leg.txt = c("Non-Zero Density", "Zero Proportion")
  leg.col = c("blue", "red")
  legend("topright", legend = leg.txt, text.col = leg.col, bty = "n")

}


