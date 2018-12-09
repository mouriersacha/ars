library("numDeriv")
library('distr')
library('testthat')
library('stats')
library('assertthat')

set.seed(1)

# Derivation function

h_x <- function(myfun) {
  fun <- function(x) log(myfun(x))
  return(fun)
}

h_x_d <- function(myfun) {
  fun <- function(x) grad(h_x(myfun), x)
  return(fun)
}

h_x_d_d <- function(myfun) {
  fun <- function(x) grad(h_x_d(myfun), x)
  return(fun)
}

# Discriminate the log_concave function

logcon <- function(myfun,lower,upper) {
  val=seq(max(lower, -10),min(upper, 10),0.1)
  val_dd=h_x_d_d(myfun)(val)
  val_dd_neg=val_dd[val_dd<0]
  assert_that(are_equal(length(val_dd), length(val_dd_neg)))
}


# Initilize the breakpoints and section

initialize <- function(myfun, lower, upper) {
  zero_point = uniroot(h_x_d(myfun),c(max(lower, -10),min(upper, 10)))$root
  if (lower==-Inf) {
    low=zero_point-1
  } else {low=(lower+zero_point)/2}
  if (upper==Inf) {
    up=zero_point+1
  } else {up=(upper+zero_point)/2}
  breakpoints = c(low,up)
  n = length(breakpoints)
  z_k = (h_x(myfun)(breakpoints[-1]) - h_x(myfun)(breakpoints[-n]) -
           breakpoints[-1] * h_x_d(myfun)(breakpoints[-1]) +
           breakpoints[-n] * h_x_d(myfun)(breakpoints[-n]))/
    (h_x_d(myfun)(breakpoints[-n]) - h_x_d(myfun)(breakpoints[-1]))
  z = c(lower,z_k,upper)
  list_return <- list(breakpoints, z)
  return(list_return)
}


# Helper Functions
# uk
uk = function(myfun,breakpoints,j) {
  uj <- function(x) h_x(myfun)(breakpoints[j]) + (x - breakpoints[j])*h_x_d(myfun)(breakpoints[j])
  return(uj)
}

# sk
sk = function(myfun,breakpoints,j) {
  sj <- function(x) exp(h_x(myfun)(breakpoints[j]) + (x - breakpoints[j])*h_x_d(myfun)(breakpoints[j]))
  return(sj)
}

# lk
lk = function(myfun,breakpoints,j) {
  n=length(breakpoints)
  if(j < n & j > 1){
    lj <- function(x) ((breakpoints[j+1] - x)*h_x(myfun)(breakpoints[j]) + (x - breakpoints[j])*h_x(myfun)(breakpoints[j+1]))/
      (breakpoints[j+1] - breakpoints[j])
  } else {lj <- function(x) -Inf }
  return(lj)
}

# z_k update

z_update <- function(myfun,breakpoints, id){
  (h_x(myfun)(breakpoints[id]) - h_x(myfun)(breakpoints[id-1]) -
     breakpoints[id] * h_x_d(myfun)(breakpoints[id]) +
     breakpoints[id-1] * h_x_d(myfun)(breakpoints[id-1]))/
    (h_x_d(myfun)(breakpoints[id-1]) - h_x_d(myfun)(breakpoints[id]))}


# Add z and breakpoints

add <- function(myfun, breakpoints, z, candidate, lower, upper) {
  id = sum(candidate > breakpoints) + 1
  breakpoints = sort(c(breakpoints,candidate))
  n = length(breakpoints)
  if(id > 0 & id < n) {
    z_new = c(z_update(myfun, breakpoints, id), z_update(myfun, breakpoints, id+1))
  } else if(id == 0){
    z_new = c(lower,z_update(myfun, breakpoints,id))
  } else if(id == n){
    z_new = c(z_update(myfun, breakpoints,id),upper)
  }

  z = sort(c(z[1:(id - 1)],z_new,z[(id+1):(n+1)]))
  return(list(breakpoints,z))
}


# sample ~ sk, first choose a sub-distribution and sample from it

sample_one <- function(myfun,breakpoints,z) {
  normalize_value = vector()
  n=length(breakpoints)
  for(i in 1:n){
    normalize_value[i] = integrate(f = sk(myfun,breakpoints,i),z[i],z[i+1])$value # when it comes to z, all changed by + 1
  }
  group = sample(1:n,size = 1,prob = normalize_value/sum(normalize_value))
  subdist <-AbscontDistribution(d=sk(myfun,breakpoints,group),low1 = max(z[group],-10000),up1 = min(z[group+1],10000))
  candidate <- r(subdist)(1)
  mylist <- list(group, candidate)
  return(mylist)
}


# Rejection

sample_one_reject <- function(myfun, breakpoints, z, samples){
  group=sample_one(myfun, breakpoints,z)[[1]]
  candidate=sample_one(myfun, breakpoints,z)[[2]]
  # w ~ U(0,1)
  w = runif(1)
  # squeezing test and rejection test (flag)

  lk_group = sum(candidate > breakpoints)
  if(w <= exp(lk(myfun,breakpoints,lk_group)(candidate) - uk(myfun,breakpoints,group)(candidate))){
    samples=c(candidate,samples)
    return(list(FALSE, samples,candidate))
  }
  else if (w <= exp(h_x(myfun)(candidate)- uk(myfun,breakpoints,group)(candidate))){
    samples=c(candidate,samples)
    return(list(TRUE, samples, candidate))
  }
  return(list(TRUE,samples,candidate))
}
