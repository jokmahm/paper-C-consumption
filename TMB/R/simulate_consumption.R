simulateConsumption <- function(sim, temp, B, par, tmin,tmax, simtemp = FALSE){
  ifelse(simtemp, 
    t <- runif(length(B), min = tmin, max = tmax),
    t <- temp)
  a <- par[1] * (t-tmin)
  h <- par[2] 
  C <- (a * B )/(1+a*h*B) + rnorm(length(B), 0, sd = par[3])
  C[C<0] <- 0 
  return(C)
}
