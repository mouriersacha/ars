library("numDeriv")
library('distr')
library('testthat')
library('stats')
library('assertthat')
set.seed(1)


# Main function
ars <- function(myfun,m,lower,upper){
  # Output
  samples = vector()
  # Whether the input function is log-concave
  logcon(myfun, lower, upper)
  # Initialize
  mylist=initialize(myfun, lower, upper)
  breakpoints=mylist[[1]]
  z=mylist[[2]]
  # Loop for samples
  Rej_test = 0
  while (length(samples) < m) {
    # update
    # print(Rej_test)
    if(Rej_test){
      mylist=add(myfun, breakpoints, z, candidate, lower, upper)
      breakpoints=mylist[[1]]
      z=mylist[[2]]
      #print(breakpoints)
    }

    # sample from sk

    mylist=sample_one_reject(myfun, breakpoints, z, samples)
    Rej_test=mylist[[1]]
    samples=mylist[[2]]
    candidate=mylist[[3]]
    if (length(samples)%%100==0){
      print(length(samples))
    }
  }
  print(length(breakpoints))
  return(samples)
}
