library("numDeriv")
library('distr')
set.seed(1)

# Main function
ars <- function(myfun,m){
  # Output
  samples = vector()
  
  # Initialize
  h_x <- function(x) log(myfun(x))
  h_x_d <- function(x) grad(h_x, x)
  
  zero_point = uniroot(h_x_d,c(-5,5))$root
  breakpoints = c(zero_point - 1, zero_point + 1)
  
  n = length(breakpoints)
  z_k = (h_x(breakpoints[-1]) - h_x(breakpoints[-n]) - 
           breakpoints[-1] * h_x_d(breakpoints[-1]) + 
           breakpoints[-n] * h_x_d(breakpoints[-n]))/
    (h_x_d(breakpoints[-n]) - h_x_d(breakpoints[-1]))
  z = c(-3,z_k,3) # boundary issue?
  
  # Helper Functions
  # uk
  uk = function(j) {
    uj <- function(x) h_x(breakpoints[j]) + (x - breakpoints[j])*h_x_d(breakpoints[j])
    return(uj)
  } 
  
  # sk
  sk = function(j) {
    sj <- function(x) exp(h_x(breakpoints[j]) + (x - breakpoints[j])*h_x_d(breakpoints[j]))
    return(sj)
  } 
  
  # lk
  lk = function(j) {
    if(j < n & j > 1){
      lj <- function(x) ((breakpoints[j+1] - x)*h_x(breakpoints[j]) + (x - breakpoints[j])*h_x(breakpoints[j+1]))/
        (breakpoints[j+1] - breakpoints[j])
    } else {lj <- function(x) -Inf }
    return(lj)
  } 
  
  # z_k update
  z_update <- function(id){
    (h_x(breakpoints[id]) - h_x(breakpoints[id-1]) - 
       breakpoints[id] * h_x_d(breakpoints[id]) + 
       breakpoints[id-1] * h_x_d(breakpoints[id-1]))/
      (h_x_d(breakpoints[id-1]) - h_x_d(breakpoints[id]))}
  
  # Loop for samples
  Rej_test = 0
  
  while (length(samples) < m) {
    # update
    if(Rej_test){
      id = sum(candidate > breakpoints) + 1
      breakpoints = sort(c(breakpoints,candidate))
      n = length(breakpoints)
      
      if(id > 0 & id < n) {
        z_new = z_update(c(id,id+1))
      } else if(id == 0){
        z_new = c(-3,z_update(id))
      } else if(id == n){
        z_new = c(z_update(id),3)
      }
      
      z = sort(c(z[1:(id - 1)],z_new,z[(id+1):(n+1)]))
    }
    
    # sample ~ sk, first choose a sub-distribution and sample from it
    normalize_value = vector()
    for(i in 1:n){
      normalize_value[i] = integrate(f = sk(i),z[i],z[i+1])$value # when it comes to z, all changed by + 1
    }
    
    group = sample(1:n,size = 1,prob = normalize_value/sum(normalize_value))
    subdist <-AbscontDistribution(d=sk(group),low1 = z[group],up1 = z[group+1])
    candidate <- r(subdist)(1)
    
    # w ~ U(0,1)
    w = runif(1)
    # squeezing test and rejection test (flag)
    lk_group = sum(candidate > breakpoints)
    Rej_test = 0
    if(w <= exp(lk(lk_group)(candidate) - uk(group)(candidate))){
      samples = c(samples,candidate)
    } else if (w <= exp(h_x(candidate)- uk(group)(candidate))){
      Rej_test = 1
      samples = c(samples,candidate)
    }
  }
  return(samples)
}