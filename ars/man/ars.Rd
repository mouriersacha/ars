\name{ars}
\alias{ars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adaptive Rejection Sampling
%%  ~~function to do ... ~~
}
\description{Method for rejection sampling for any univariate log-concave probability density function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ars(myfun, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myfun}{density the user wants to sample from
%%     ~~Describe \code{x} here~~
}
  \item{m}{sample size}
}
\details{The method is adaptive: as sampling proceeds, the rejection envelope and the squeezing function converge to the density function. The rejection envelope and squeezing function are piecewise exponential functions forming arcs between those points of contact. The technique is intended for situations where evaluation of the density is computationally expensive.
%%  ~~ If necessary, more details than the description above ~~
}
\value{a sample of inputted length n from a density
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Gilks, W. R, Wild, P. (1992) Adaptive Rejection Sampling for Gibbs Sampling, \emph{Applied Statistics} 41:337-348
%% ~put references to the literature/web site here ~
}
\author{Cindy Zhang, Sihan Chen, Sacha Mourier
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(ars)

#Example 1: sample 500 values from the standard normal distribution
ars1 <- ars(myfun = dnorm, m = 500)
hist(ars1, breaks = 100)

#Example 2: sample 400 values from the t-distribution
t_dist <- function (x) dt(x, df = 12)
ars2 <- ars(myfun = t_distr, m = 400)
hist(ars2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
