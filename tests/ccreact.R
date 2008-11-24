if(!require("BB"))stop("this test requires package BB.")
iseed <- 1236  # this seed was used for tests conducted on March 25, 2008.  
set.seed(iseed)

ccreact <- function(x) {
n <- length(x)
f <- rep(NA, n)
a <- sqrt(2) - 1
j <- 4:(n-2)
f[1] <- x[1] - (1 - x[1])* x[3] - a * (1 + 4*x[2])
f[2] <- - (1 - x[1])* x[4] - a * (1 + 4*x[2])
f[3] <- a * x[1] - (1 - x[1])* x[5] - x[3] * (1 + 4*x[4])
f[j] <- x[1] * x[j-2] - (1 - x[1]) * x[j+2] - x[j] * (1 + 4*x[j-1])
f[n-1] <- x[1] * x[n-3] - x[n-1] * (1 + 4*x[n-2])
f[n] <- x[1] * x[n-2] - (1 - x[1]) - x[n] * (1 + 4*x[n-1])
f
}


p0 <- rep(c(0.1,0.2,0.3,0.4,0.5,0.4,0.3,0.2), 5)
#ans <- nlsolve(p0, fn=ccreact, method="L-BFGS-B")

ans1 <- dfsane(par=p0, fn=ccreact, method=1)
ans2 <- dfsane(par=p0, fn=ccreact, method=2)
ans3 <- sane(par=p0, fn=ccreact, method=1)
ans4 <- sane(par=p0, fn=ccreact, method=2)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


p1 <- sort(runif(5))
p0 <- rep(c(p1, rev(p1)[-c(1,5)]), 5)
#nlsolve(p0, ccreact)

ans1 <- dfsane(par=p0, fn=ccreact, method=1)
ans2 <- dfsane(par=p0, fn=ccreact, method=2, control = list(NM=TRUE))
ans3 <- sane(par=p0, fn=ccreact, method=1)
ans4 <- sane(par=p0, fn=ccreact, method=2)

z <- c(ans1$resid, ans2$resid, ans3$resid, ans4$resid)
print(z, digits=18)
if(any(abs(
 c(0.1676742849181734, 0.0890678685479564, 0.1676742849181734, 0.1699072249956761)
    - z) > 1e-10)) stop("test ccreact xx FAILED")

c(ans1$feval, ans2$feval, ans3$feval, ans4$feval) 

