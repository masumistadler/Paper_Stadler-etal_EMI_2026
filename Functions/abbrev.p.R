# Significance function
abbrev.p <- function(x){
  if(x < 0.0001){
    out <- c("< 0.0001","****")
  } else if(x <= 0.001){
    out <- c("< 0.001","***")
  } else if(x <= 0.01){
    out <- c("< 0.01","**")
  } else if(x <= 0.05){
    out <- c("< 0.05","*")
  } else {
    out <- c(as.character(paste("=",round(x, 2))),"n.s.")
  }
  return(out)
}