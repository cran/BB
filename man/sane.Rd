\name{sane}
\alias{sane}
\title{ Solving Large-Scale Nonlinear System of Equations}
\description{Non-Monotone spectral approach for Solving Large-Scale Nonlinear Systems of Equations}
\usage{
  sane(par, fn, method=1, control=list(), ...) 
 }
\arguments{
    \item{fn}{a function that takes a real vector as argument and  returns a real vector of same length (see details).}
 
    \item{par}{A real vector argument to \code{fn}, indicating the 
      initial guess for the root of the nonlinear system.}

     \item{method}{An integer (1, 2, or 3) specifying which Barzilai-Borwein steplength to use.  The default is 1.  See *Details*.} 

    \item{control}{A list of control parameters. See *Details*.}
 
    \item{...}{Additional arguments passed to \code{fn}.}
}
\value{ A list with the following components:
\item{par}{The best set of parameters that solves the nonlinear system.}

\item{residual}{L2-norm of the function evaluated at \code{par}, divided
  by \code{sqrt(npar)}, where \code{npar} is the number of parameters.}

\item{fn.reduction}{Reduction in the L2-norm of the function from the initial
  L2-norm.}

\item{feval}{Number of times \code{fn} was evaluated.}

\item{iter}{Number of iterations taken by the algorithm.}

\item{convergence}{An integer code indicating type of convergence.  \code{0} 
  indicates successful convergence, in which case the \code{resid} is smaller 
  than \code{tol} Error codes are \code{1} indicates that the iteration 
  limit \code{maxit} has been reached.  \code{2} indicates failure due to an 
  error in function evaluation.  \code{3} indicates failure due to 
  exceeding 100 steplength reductions in line-search.  \code{4} indicates 
  failure due to an anomalous iteration.}

\item{message}{A text message explaining which termination criterion was used.}
}

\details{
 The function \code{sane} implements a non-monotone spectral residual method for finding a root of nonlinear systems.  It stands for "spectral approach for nonlinear equations".  
It differs from the function \code{dfsane} in that it requires an approximation of a directional derivative at every iteration of the merit function $F(x)^t F(x)$.

R adaptation, with significant modifications, by Ravi Varadhan, Johns Hopkins 
University (March 25, 2008), from the original FORTRAN code of La Cruz and 
Raydan (2003) 
available at \url{http://kuainasi.ciens.ucv.ve/ccct/mraydan/mraydan.html}.  
.  

A major modification in our R adaptation of the original FORTRAN code is the availability of 3 different options for Barzilai-Borwein (BB) steplengths: \code{method = 1} is the BB 
steplength used in LaCruz and Raydan (2003);  \code{method = 2} is equivalent to the other steplength proposed in Barzilai and Borwein's (1988) original paper.  
Finally, \code{method = 3}, is a new steplength, which is equivalent to that first proposed in Varadhan and Roland (2008) for accelerating the EM algorithm.  
In fact, Varadhan and Roland (2008) considered 3 equivalent steplength schemes in their EM acceleration work.  Here, we have chosen \code{method = 1} 
as the "default" method.  However, we have not seen major differences between the three steplength schemes in our experiments.  

Argument \code{control} is a list specifing any changes to default values of 
algorithm control parameters.  Note that the names of these must be 
specified completely.  Partial matching will not work.

\item{M}{A positive integer, typically between 5-20, that controls the monotonicity of the algorithm.  \code{M=1} would enforce strict monotonicity 
in the reduction of L2-norm of \code{fn}, whereas larger values allow for more non-monotonicity.  Global convergence under non-monotonicity is ensured by 
enforcing the Grippo-Lampariello-Lucidi condition (Grippo et al. 1986) in a non-monotone line-search algorithm.  Values of \code{M} between 5 to 20 are generally good.  The default is \code{M = 10}.}  

\item{maxit}{The maximum number of iterations.  The default is \code{maxit = 1500}.}

\item{tol}{The absolute convergence tolerance on the residual L2-norm 
  of \eqn{F}. Convergence is declared  
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

   W LaCruz, and M Raydan (2003),  Nonmonotone spectral methods for large-scale 
   nonlinear systems, \emph{Optimization Methods and Software}, 18, 583-599
   (see \url{http://kuainasi.ciens.ucv.ve/ccct/mraydan/mraydan.html}). 

   R Varadhan and C Roland (2008),  Simple and globally-convergent methods for accelerating the convergence of any EM algorithm,  \emph{Scandinavian J Statistics}, doi: 10.1111/j.1467-9469.2007.00585.x.  

   R Varadhan and PD Gilbert (2008),  BB: An R package of Barzilai-Borwein spectral methods for solving and optimizing large-scale nonlinear systems, Unpublished.

}

\seealso{
  \code{\link{dfsane}},
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
    sane(par=p0, fn=trigexp)
    sane(par=p0, fn=trigexp, method=2)
    sane(par=p0, fn=trigexp, control=list(triter=5, M=20))
}

\keyword{multivariate}
