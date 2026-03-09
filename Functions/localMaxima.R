# Identify positions of local maxima
# this function is very sensitive some quality control has to be done afterwards
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(-Inf, x)) > 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples the every second location (i.e. switch to TRUE)
  
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}