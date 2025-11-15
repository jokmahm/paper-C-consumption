simulateConsumption <- function(sim, temp, B, P, par, tmin,tmax, simtemp = FALSE){
  ifelse(simtemp, 
    t <- runif(length(B), min = tmin, max = tmax),
    t <- temp)
  a <- par[2] * (t-tmin)*(tmax-t)^(par[4])
  h <- par[3] 
  C <- (a * B*P^(-par[1]) )/(1+a*h*B*P^(-par[1])) + rnorm(length(B), 0, sd = par[5])
  C[C<0] <- 0 
  return(C)
}
