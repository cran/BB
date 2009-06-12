\name{multiStart}
\alias{multiStart}
\title{Nonlinear Optimization or Root-Finding with Multiple Starting Values}
\description{Start \code{BBsolve} or \code{BBoptim} from multiple starting 
points to obtain multiple solutions and to test sensitivity to starting values.}
\usage{
  multiStart(par, fn, gr=NULL, action = c("solve", "optimize"), 
	method=c(2,3,1), control=list(), details=FALSE, ...) 
  }
\arguments{
    \item{par}{A real matrix, each row of which is an argument to \code{fn},
      indicating initial guesses for solving a nonlinear system \code{fn = 0} 
      or for optimizing the objective function \code{fn}.}
    \item{fn}{see \code{BBsolve} or \code{BBoptim}.} 
    \item{gr}{Only required for optimization.  See \code{BBoptim}.} 
    \item{action}{A character string indicating whether to solve a nonlinear 
      system or to optimize.  Default is ``solve''.} 
    \item{method}{see \code{BBsolve} or \code{BBoptim}.} 
    \item{control}{See \code{BBsolve} and \code{BBoptim}.}
    \item{details}{Logical indicating if the result should include the full
      result from \code{BBsolve} or \code{BBoptim} for each starting value.}
    \item{...}{arguments passed fn (via the optimization algorithm).}
}
\value{list with elements \code{par}, \code{values}, and \code{converged}. 
  It optionally returns an attribute called ``details'', which is a list as long as 
  the number of starting values, which contains the complete object returned 
  by \code{dfsane} or \code{spg} for each starting value.}

\details{
The optimization or root-finder is run with each row of \code{par} indicating 
initial guesses.
}
\seealso{
  \code{\link{BBsolve}},
  \code{\link{BBoptim}},
  \code{\link{dfsane}},
  \code{\link{spg}}
   }
   
\examples{
# Use a preset seed so the example is reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=1234))

# Finding multiple roots of a nonlinear system
brownlin <- function(x) {
# Brown's almost linear system(A.P. Morgan, ACM 1983)
# two distinct solutions if n is even
# three distinct solutions if n is odd  
  	n <- length(x)
  	f <- rep(NA, n)
	nm1 <- 1:(n-1)
	f[nm1] <- x[nm1] + sum(x) - (n+1)
	f[n] <- prod(x) - 1 
	f
}

p <- 9
n <- 100
p0 <- matrix(rnorm(n*p), n, p)  # n starting values, each of length p
ans <- multiStart(par=p0, fn=brownlin)
pmat <- ans$par[ans$conv, 1:p] # selecting only converged solutions
ord1 <- order(abs(pmat[,1]))
round(pmat[ord1, ], 3)  # all 3 roots can be seen

# An optimization example
rosbkext <- function(x){
n <- length(x)
j <- 2 * (1:(n/2))
jm1 <- j - 1
sum(100 * (x[j] - x[jm1]^2)^2 + (1 - x[jm1])^2)
}

p0 <- rnorm(50)
spg(par=p0, fn=rosbkext)
BBoptim(par=p0, fn=rosbkext)

pmat <- matrix(rnorm(500), 100, 5)  # 100 starting values each of length 5 
ans <- multiStart(par=pmat, fn=rosbkext, action="optimize")
ans
attr(ans, "details")[[1]]  # 

pmat <- ans$par[ans$conv, 1:5] # selecting only converged solutions
round(pmat, 3)
}

\keyword{multivariate}
