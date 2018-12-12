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

# Discriminate the positive function

pos <- function(myfun,lower,upper) {
  val=seq(max(lower, -10),min(upper, 10),0.1)
  val=myfun(val)
  val_pos=val[val>0]
  assert_that(are_equal(length(val), length(val_pos)))
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
  breakpoints_d=c(h_x_d(myfun)(breakpoints[1]), h_x_d(myfun)(breakpoints[2]))
  list_return <- list(breakpoints, breakpoints_d, z)
  return(list_return)
}

initialize(dnorm, -Inf, 5)


# Helper Functions
# uk
uk = function(myfun,breakpoints,breakpoints_f,breakpoints_d,j) {
  uj <- function(x) breakpoints_f[j] + (x - breakpoints[j])*breakpoints_d[j]
  return(uj)
}

# sk
sk = function(myfun,breakpoints,breakpoints_f,breakpoints_d,j) {
  sj <- function(x) exp(breakpoints_f[j] + (x - breakpoints[j])*breakpoints_d[j])
  return(sj)
}

# lk
lk = function(myfun,breakpoints,breakpoints_f,j) {
  n=length(breakpoints)
  if(j < n & j >= 1){
    lj <- function(x) ((breakpoints[j+1] - x)*breakpoints_f[j] + (x - breakpoints[j])*breakpoints_f[j+1])/
      (breakpoints[j+1] - breakpoints[j])
  } else {lj <- function(x) -Inf }
  return(lj)
}

# z_k update

z_update <- function(myfun,breakpoints,breakpoints_f, breakpoints_d, id){
  (breakpoints_f[id] - breakpoints_f[id-1] -
     breakpoints[id] * breakpoints_d[id] +
     breakpoints[id-1] * breakpoints_d[id-1])/
    (breakpoints_d[id-1] - breakpoints_d[id])}


# Add z and breakpoints

add <- function(myfun, breakpoints,breakpoints_f, breakpoints_d, z, candidate, lower, upper) {
  id = sum(candidate > breakpoints) + 1
  n=length(breakpoints)
  if (id!=(n+1)) {
    for (i in (n+1):(id+1)) {
      breakpoints[i]=breakpoints[i-1]
      breakpoints_f[i]=breakpoints_f[i-1]
      breakpoints_d[i]=breakpoints_d[i-1]
    }
  }
  #print(breakpoints)
  breakpoints[id]=candidate
  breakpoints_f[id]=h_x(myfun)(candidate)
  breakpoints_d[id]=h_x_d(myfun)(candidate)
  n = n+1
  if(id > 1 & id < n) {
    z_new = c(z_update(myfun, breakpoints, breakpoints_f, breakpoints_d, id), z_update(myfun,breakpoints, breakpoints_f, breakpoints_d, id+1))
    z = c(z[1:(id - 1)],z_new,z[(id+1):n])
  } else if(id == 1){
    z_new = c(lower,z_update(myfun, breakpoints, breakpoints_f, breakpoints_d,id+1))
    z = c(z_new,z[(id+1):n])
  } else if(id == n){
    z_new = c(z_update(myfun, breakpoints, breakpoints_f, breakpoints_d, id),upper)
    z = c(z[1:(id - 1)],z_new)
  }
  return(list(breakpoints,breakpoints_f,breakpoints_d,z,id))
}


# sample ~ sk, first choose a sub-distribution and sample from it

sample_one <- function(myfun,normalize_value, breakpoints,breakpoints_f,breakpoints_d,z) {
  n=length(breakpoints)
  group = sample(1:n,size = 1,prob = normalize_value/sum(normalize_value))
  subdist <-AbscontDistribution(d=sk(myfun,breakpoints,breakpoints_f,breakpoints_d,group),low1 = max(z[group],-10000),up1 = min(z[group+1],10000))
  candidate <- r(subdist)(1)
  mylist <- list(group, candidate)
  return(mylist)
}


# Rejection

sample_one_reject <- function(myfun, normalize_value, breakpoints, breakpoints_f,breakpoints_d,z, samples){
  mylist=sample_one(myfun, normalize_value, breakpoints,breakpoints_f,breakpoints_d,z)
  group=mylist[[1]]
  candidate=mylist[[2]]
  # w ~ U(0,1)
  w = runif(1)
  # squeezing test and rejection test (flag)

  lk_group = sum(candidate > breakpoints)
  if(w <= exp(lk(myfun,breakpoints,breakpoints_f,lk_group)(candidate) - uk(myfun,breakpoints,breakpoints_f,breakpoints_d,group)(candidate))){
    samples=c(candidate,samples)
    return(list(FALSE, samples,candidate))
  }
  else if (w <= exp(h_x(myfun)(candidate)- uk(myfun,breakpoints,breakpoints_f,breakpoints_d,group)(candidate))){
    samples=c(candidate,samples)
    return(list(TRUE, samples, candidate))
  }
  return(list(TRUE,samples,candidate))
}

