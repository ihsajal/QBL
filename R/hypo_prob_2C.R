Sprob <- function(z,nz,a,b,e){
  prob <- 0
  prob <- (((b+e)^a - b^a)/((b+e)^a))^z
  return(prob)
}

