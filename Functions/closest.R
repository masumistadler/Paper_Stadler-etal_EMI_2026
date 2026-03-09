## Finds the above and below closest value
closest <- function(x,value){
  b <- x[x < value]
  if(length(b) == 0){
    below <- NA
  } else {
    below <- b[which.min(abs(b-value))]
  }
  
  a <- x[x > value]
  if(length(a) == 0){
    above <- NA
  } else {
    above <- a[which.min(abs(a-value))]
  }
  
  if(length(c(a,b)) == 0L){
    above <- x; below <- x
  }
  
  return(c(below, above))
}