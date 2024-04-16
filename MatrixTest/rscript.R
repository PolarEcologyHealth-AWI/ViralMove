Mat4Array      <- function(x) {
  
  rows = length(x)
  cols = length(x[[1]])
  d3 = length(x[[1]][[1]])
  d4 = length(x[[1]][[1]][[1]])
  
  out <- array(dim = c(rows, cols, d3, d4))
  
  for(r in 1:rows) {
    for(c in 1:cols) {
      for(d.1 in 1:d3) {
        for(d.2 in 1:d4) {
        out[r, c, d.1, d.2] <- x[[r]][[c]][[d.1]][[d.2]]
        }
      }
    }
  }
  out 
}


Mat3Array(timesTwo(2,3,4,2))