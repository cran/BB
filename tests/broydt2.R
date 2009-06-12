if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1236)
old.seed <- setRNG(test.rng)
#iseed <- 1236  # this seed was used for tests conducted on March 25, 2008.  
#set.seed(iseed)

broydt <- function(x) {
n <- length(x)
f <- rep(NA, n)
#h <- 0.5
h <- 2
f[1] <- ((3 - h*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - h*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - h*x[n]) * x[n]) - x[n-1] + 1
f
}

p0 <- -runif(500)
ans1 <- dfsane(par=p0, fn=broydt, method=1)
ans2 <- dfsane(par=p0, fn=broydt, method=2)
ans3 <- sane(par=p0, fn=broydt, method=2)
ans4 <- sane(par=p0, fn=broydt, method=3)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 

nsim <- 100
dfsane1.broydt <- dfsane2.broydt <- sane1.broydt <- sane2.broydt <- matrix(NA, nsim, 5)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- -runif(500)
t1 <- system.time(ans <- sane(par=p0, fn=broydt, method=1, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))sane1.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t1)
t2 <- system.time(ans <- sane(par=p0, fn=broydt, method=2, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))sane2.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t2)
t3 <- system.time(ans <- dfsane(par=p0, fn=broydt, method=1, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))dfsane1.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t3)
t4 <- system.time(ans <- dfsane(par=p0, fn=broydt, method=2, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans)) dfsane2.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t4)
}


z <- apply(sane1.broydt, 2, summary)
print(z)
print(z[,1], digits=18)
if(any(abs(c(5.278e-08, 6.382e-08, 6.969e-08, 7.106e-08, 7.653e-08, 9.932e-08)
    - z[,1]) > 1e-10)) stop("test broydt2 sane1.broydt FAILED")

z <- apply(sane2.broydt, 2, summary)
print(z)
print(z[,1], digits=18)
if(any(abs(c(5.612e-08, 6.538e-08, 7.172e-08, 7.267e-08, 8.026e-08, 9.912e-08) 
    - z[,1]) > 1e-10)) stop("test broydt2 sane2.broydt FAILED")

z <- apply(dfsane1.broydt, 2, summary)
print(z)
print(z[,1], digits=18)
if(any(abs(c(5.278e-08, 6.382e-08, 6.969e-08, 7.106e-08, 7.653e-08, 9.932e-08)  
    - z[,1]) > 1e-10)) stop("test broydt2 dfsane1.broydt FAILED")

z <- apply(dfsane2.broydt, 2, summary)
print(z)
print(z[,1], digits=18)
if(any(abs(c(5.612e-08, 6.538e-08, 7.172e-08, 7.267e-08, 8.026e-08, 9.912e-08)  
    - z[,1]) > 1e-10)) stop("test broydt2 dfsane2.broydt FAILED")

