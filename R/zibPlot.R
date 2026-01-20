#' Visualise a zero-inflated beta regression model for simulated data
#'
#' Plots a simulated zero-inflated beta regression dataset of class
#' \code{"zibSim"}, separately displaying the zero abundance in the dataset and
#' fitting a beta regression model to the distribution of the non-zero data,
#' which is displayed on a histogram.
#'
#' @param x An object of class \code{"zibSim"}.
#' @param microbe_name Name of simulated microbe to title the plot; "Microbe" by
#' default.
#' @param ... Catches unused arguments to \code{plot} (not currently
#' implemented).
#'
#' @export
#' @importFrom MASS "fitdistr"
#'
#' @seealso \code{\link{plot}}, \code{\link{zibSim}}
#'
#' @examples
#' n = 1000
#' t = 3
#' X <- as.matrix(c(rep(0, n/2 * t), rep(1, n/2 * t)))
#'
#' y1 <- zibSim(n = n, t = t, a = 0.5, b = 0.8, sigma1 = 2.5, sigma2 = 0.5,
#' phi = 7.2, alpha = 0.5, beta = -0.5, X = X, Z = X, seed = 6874)
#'
#' plot(y1, microbe_name = "Simulated Data")
#'

plot.zibSim <- function(x, ..., microbe_name = "Microbe") {

  # get the data
  df <- x$rel_abundance

  # create the plot
  # par(new = TRUE)
  breaks <- seq(0, 1, by = 0.025)
  non_zero <- subset(df, df != 0)

  fit <- MASS::fitdistr(non_zero, stats::dbeta, start = list(shape1 = 1, shape2 = 1))

  zero_prop <- mean(df == 0)
  graphics::par(mar = c(5, 4, 4, 5) + 0.1)
  graphics::par(col.axis = "blue", col.lab = "blue")
  graphics::hist(df, freq = FALSE, xlim = c(0, 1), breaks = breaks,
       main = microbe_name, xlab = "Abundance", ylab = "Density")
  graphics::curve(stats::dbeta(x, fit$estimate["shape1"], fit$estimate["shape2"]),
        col = "blue", lwd = 2, add = TRUE)
  graphics::par(new = TRUE, col.axis = "red")
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = range(df), ylim = c(0, 1))
  graphics::lines(c(0, 0), c(0, zero_prop), col = "red", lwd = 3)
  graphics::axis(side = 4, at = seq(0, 1, 0.2))
  graphics::mtext("Zero proportion", side = 4, line = 3, col = "red")
  leg.txt = c("Non-Zero Density", "Zero Proportion")
  leg.col = c("blue", "red")
  graphics::legend("topright", legend = leg.txt, text.col = leg.col, bty = "n")

}
