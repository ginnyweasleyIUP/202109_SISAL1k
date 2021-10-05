LPlot_try <- function (x, conf = TRUE, col = "black", alpha = 0.3, removeFirst = 0, 
          removeLast = 0, xlab = "f", ylab = "PSD", ...) 
{
  index <- (removeFirst + 1):(length(x$freq) - removeLast)
  x$freq <- x$freq[index]
  x$spec <- x$spec[index]
  x$lim.1 <- x$lim.1[index]
  x$lim.2 <- x$lim.2[index]
  plot(x$freq, x$spec, type = "l", col = col, xlab = xlab, 
       ylab = ylab, ...)
  if (conf) 
    polygon(c(x$freq, rev(x$freq)), c(x$lim.1, rev(x$lim.2)), 
            col = ColTransparent(col, alpha), border = NA)
  lines(x$freq, x$spec, col = col, ...)
}