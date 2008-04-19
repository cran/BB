\name{dfsane}
\alias{dfsane}
\title{ Solving Large-Scale Nonlinear System of Equations}
\description{Derivative-Free Spectral Approach for solving nonlinear systems of equations}
\usage{
  dfsane(par, fn, method=2, control=list(), ...) 
 }
\arguments{
    \item{fn}{a function that takes a real vector as argument and  returns a real vector of same length (see details).}
 
    \item{par}{A real vector argument to \code{fn}, indicating the 
      initial guess for the root of the nonlinear system.}
 
    \item{method}{An integer (1, 2, or 3) specifying which Barzilai-Borwein steplength to use.  The default is 2.  See *Details*.} 
 
    \item{control}{A list of control parameters. See *Details*.}
 
    \item{...}{Additional arguments passed to \code{fn}.} 
}
\value{A list with the following components:
\item{par}{The best set of parameters that solves the nonlinear system.}

\item{residual}{L2-norm of the function at convergence, divided by sqrt(npar),
  where "npar" is the number of parameters.}

\item{fn.reduction}{Reduction in the L2-norm of the function from the initial
  L2-norm.}

\item{feval}{Number of times \code{fn} was evaluated.}

\item{iter}{Number of iterations taken by the algorithm.}

\item{convergence}{An integer code indicating type of convergence.  \code{0} 
  indicates successful convergence, in which case the \code{resid} is smaller 
  than \code{tol}. Error codes are \code{1} indicates that the iteration 
  limit \code{maxit} has been reached.  \code{2} indicates failure due to an 
  error in function evaluation.  \code{3} indicates failure due to 
  exceeding 100 steplength reductions in line-search.}  

\item{message}{A text message explaining which termination criterion was used.}
}

\details{
The function \code{dfsane} is another algorithm for implementing non-monotone spectral residual method for finding a root of nonlinear systems, by working without gradient information.  
It stands for "derivative-free spectral approach for nonlinear equations". It differs from the function \code{sane} in that \code{sane} requires an approximation 
of a directional derivative at every iteration of the merit 
function \eqn{F(x)^t F(x)}, where \eqn{F} is specified by the argument \code{fn}.

R adaptation, with significant modifications, by Ravi Varadhan, Johns Hopkins University (March 25, 2008), from the original FORTRAN code of La Cruz, Martinez, and Raydan (2006).  

A major modification in our R adaptation of the original FORTRAN code is the availability of 3 different options for Barzilai-Borwein (BB) steplengths: \code{method = 1} is the BB 
steplength used in LaCruz, Martinez and Raydan (2006);  \code{method = 2} is equivalent to the other steplength proposed in Barzilai and Borwein's (1988) original paper.  
Finally, \code{method = 3}, is a new steplength, which is equivalent to that first proposed in Varadhan and Roland (2008) for accelerating the EM algorithm.  
In fact, Varadhan and Roland (2008) considered 3 similar steplength schemes in their EM acceleration work.  Here, we have chosen \code{method = 2} 
as the "default" method, since it generally performs better than others.  However, it is not universally better and in some cases, one of the other two 
steplength schemes perform better.  Therefore, if the default does not perform well in a particular  nonlinear system, we recommend trying the other two steplength schemes. 

Argument \code{control} is a list specifing any changes to default values of algorithm control parameters.  Note that the names of these must be 
specified completely.  Partial matching does not work.
The list items are as follows:

\item{M}{A positive integer, typically between 5-20, that controls the monotonicity of the algorithm.  \code{M=1} would enforce strict monotonicity 
in the reduction of L2-norm of \eqn{F}, whereas larger values allow for more non-monotonicity.  Global convergence under non-monotonicity is ensured by 
enforcing the Grippo-Lampariello-Lucidi condition (Grippo et al. 1986) in a non-monotone line-search algorithm.  Values of \code{M} 
between 5 to 20 are generally good.  The default is \code{M = 10}.} 

\item{maxit}{The maximum number of iterations.  The default is \code{maxit = 1500}.}

\item{tol}{The absolute convergence tolerance on the residual L2-norm 
  of \eqn{F}.  Convergence is declared 
  when \eqn{\|F(x)\| / \sqrt(npar) < \mbox{tol}}{sum(F(x)^2) / npar < tol^2}.  
Default is \code{tol = 1.e-07}.}

\item{trace}{A logical variable (TRUE/FALSE).  If \code{TRUE}, information on the progress of solving the system is produced.  
Default is \code{trace = TRUE}.}

\item{triter}{An integer that controls the frequency of tracing when \code{trace=TRUE}. Default is \code{triter=10}, which means that
the L2-norm of \code{fn} is printed at every 10-th iteration.}
}

\references{ 

   J Barzilai, and JM Borwein  (1988),  Two-point step size gradient methods, \emph{IMA J Numerical Analysis}, 8, 141-148. 

   L Grippo, F Lampariello, and S Lucidi (1986),  A nonmonotone line search technique for Newton's method,  \emph{SIAM J on Numerical Analysis}, 23, 707-716. 

   W LaCruz, JM Martinez, and M Raydan (2006),  Spectral residual method without gradient information for solving large-scale nonlinear systems of equations, \emph{Mathematics of Computation}, 75, 1429-1448. 

   R Varadhan and C Roland (2008),  Simple and globally-convergent methods for accelerating the convergence of any EM algorithm,  \emph{Scandinavian J Statistics}, doi: 10.1111/j.1467-9469.2007.00585.x. 

   R Varadhan and PD Gilbert (2008),  BB: An R package of Barzilai-Borwein spectral methods for solving and optimizing large-scale nonlinear systems,  \emph{Unpublished}. 

}

\seealso{
  \code{\link{sane}},
  \code{\link{spg}},
  \code{\link[numDeriv]{grad}}
   }

\examples{
  trigexp <- function(x) {
# Test function No. 12 in the Appendix of LaCruz and Raydan (2003)
    n <- length(x)
    F <- rep(NA, n)
    F[1] <- 3*x[1]^2 + 2*x[2] - 5 + sin(x[1] - x[2]) * sin(x[1] + x[2])
    tn1 <- 2:(n-1)
    F[tn1] <- -x[tn1-1] * exp(x[tn1-1] - x[tn1]) + x[tn1] * ( 4 + 3*x[tn1]^2) +
        2 * x[tn1 + 1] + sin(x[tn1] - x[tn1 + 1]) * sin(x[tn1] + x[tn1 + 1]) - 8 
    F[n] <- -x[n-1] * exp(x[n-1] - x[n]) + 4*x[n] - 3
    F
    }

    p0 <- rnorm(1000)
    dfsane(par=p0, fn=trigexp)  # default is method=2
    dfsane(par=p0, fn=trigexp, method=1)
    dfsane(par=p0, fn=trigexp, method=3)
    dfsane(par=p0, fn=trigexp, control=list(triter=5, M=5))
}
\keyword{multivariate}
