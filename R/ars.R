library("numDeriv")
library('distr')
library('testthat')
library('stats')
library('assertthat')

set.seed(1)


# Main function
ars <- function(myfun,m,lower,upper){
  n=nargs()
  assert_that(are_equal(n,4))
  assert_that(is.function(myfun))
  # Output
  samples = vector()
  # Whether the input function is log-concave
  assert_that(is.numeric(lower))
  assert_that(is.numeric(upper))
  assert_that(upper>lower)

  pos(myfun, lower, upper)
  logcon(myfun, lower, upper)
  # Initialize
  mylist=initialize(myfun, lower, upper)
  breakpoints=mylist[[1]]
  breakpoints_f=breakpoints_f=c(h_x(myfun)(breakpoints[1]), h_x(myfun)(breakpoints[2]))
  breakpoints_d=mylist[[2]]
  z=mylist[[3]]
  # normalize_value
  n=length(breakpoints)
  normalize_value=vector()
  for(i in 1:n){
    normalize_value[i] = exp(breakpoints_f[i]-breakpoints[i]*breakpoints_d[i])/breakpoints_d[i]*(exp(breakpoints_d[i]*z[i+1])-exp(breakpoints_d[i]*z[i])) # when it comes to z, all changed by + 1
  }
  # Loop for samples
  Rej_test = 0
  while (length(samples) < m) {
    # update
    # print(Rej_test)
    if(Rej_test){
      mylist=add(myfun, breakpoints,breakpoints_f,breakpoints_d, z, candidate, lower, upper)
      breakpoints=mylist[[1]]
      breakpoints_f=mylist[[2]]
      breakpoints_d=mylist[[3]]
      z=mylist[[4]]
      id=mylist[[5]]
      if (id<n) {
        for (i in (n+1):(id+2)) {
          normalize_value[i]=normalize_value[i-1]
        }
      }
      if (id!=(n+1)) {
        normalize_value[id+1]=exp(breakpoints_f[id+1]-breakpoints[id+1]*breakpoints_d[id+1])/breakpoints_d[id+1]*(exp(breakpoints_d[id+1]*z[id+2])-exp(breakpoints_d[id+1]*z[id+1]))
      }
      normalize_value[id]=exp(breakpoints_f[id]-breakpoints[id]*breakpoints_d[id])/breakpoints_d[id]*(exp(breakpoints_d[id]*z[id+1])-exp(breakpoints_d[id]*z[id]))
      if (id!=1){
        normalize_value[id-1]=exp(breakpoints_f[id-1]-breakpoints[id-1]*breakpoints_d[id-1])/breakpoints_d[id-1]*(exp(breakpoints_d[id-1]*z[id])-exp(breakpoints_d[id-1]*z[id-1]))

      }
      n=n+1
      #print(breakpoints)
    }

    # sample from sk
    #print(normalize_value)
    mylist=sample_one_reject(myfun, normalize_value, breakpoints, breakpoints_f,breakpoints_d,z, samples)
    Rej_test=mylist[[1]]
    samples=mylist[[2]]
    candidate=mylist[[3]]
  }
  #print(length(breakpoints))
  return(samples)
}


