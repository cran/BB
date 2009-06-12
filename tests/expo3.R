if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1236)
old.seed <- setRNG(test.rng)
#iseed <- 1236  # this seed was used for tests conducted on March 25, 2008.  
#set.seed(iseed)

expo3 <- function(p) {
#  From La Cruz and Raydan, Optim Methods and Software 2003, 18 (583-599)
n <- length(p)
f <- rep(NA, n)
onm1 <- 1:(n-1) 
f[onm1] <- onm1/10 * (1 - p[onm1]^2 - exp(-p[onm1]^2))
f[n] <- n/10 * (1 - exp(-p[n]^2))
f
}

n <- 200
#p0 <- (1:n)/(4*n^2)
p0 <- rnorm(n, sd=2)
ans1 <- dfsane(par=p0, fn=expo3, method=1)
ans2 <- dfsane(par=p0, fn=expo3, method=2)
ans3 <- sane(par=p0, fn=expo3)
ans4 <- sane(par=p0, fn=expo3, method=3)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 

# switched BFGS=TRUE to BFGS=FALSE below for speed, and 500 to 100
nsim <- 100
dfsane1.expo3 <- dfsane2.expo3 <- sane1.expo3 <- sane2.expo3 <- matrix(NA, nsim, 5)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(500)
t1 <- system.time(ans <- sane(par=p0, fn=expo3, method=1,
                   control=list(BFGS=FALSE, trace=FALSE)))[1]
if (!is.null(ans))sane1.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t1)
t2 <- system.time(ans <- sane(par=p0, fn=expo3, method=2,
                    control=list(BFGS=FALSE, trace=FALSE)))[1]
if (!is.null(ans))sane2.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t2)
t3 <- system.time(ans <- dfsane(par=p0, fn=expo3, method=1,
                    control=list(BFGS=FALSE, trace=FALSE)))[1]
if (!is.null(ans))dfsane1.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t3)
t4 <- system.time(ans <- dfsane(par=p0, fn=expo3, method=2,
                    control=list(BFGS=FALSE, trace=FALSE)))[1]
if (!is.null(ans)) dfsane2.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t4)
}

z <- apply(sane1.expo3, 2, summary)
print(z)
print(z[,1], digits=16)
good <- c(
    0,          0,        0,      1.294e+00, 3.898e-01, 3.189e+01)  

#5.043e-09  5.475e-08  9.066e-08  1.303e+00  3.898e-01  3.189e+01 # Linux 64
#1.445e-08  4.596e-08  9.197e-08  1.294e+00  3.898e-01  3.189e+01 # Linux 32
#9.597e-09  6.499e-08  9.594e-08  1.294e+00  3.898e-01  3.189e+01 # CRAN Win-builder
print(good - z[,1], digits=16)
if(any(abs(good - z[,1]) * c(1,1,1,1e-6,1,1) > 1))# 1e-8))
     stop("test expo3 sane1.expo3 FAILED")


z <- apply(sane2.expo3, 2, summary)
print(z)
print(z[,1], digits=16)
good <- c(
    0,          0,        0,      3.402e-02, 8.052e-08, 2.106e+00)  
#1.443e-09  1.736e-08  5.522e-08  3.402e-02  8.052e-08  2.106e+00 # Linux 32
#2.041e-09  2.685e-08  5.552e-08  3.258e-02  8.407e-08  2.105e+00 # Linux 64
print(good - z[,1], digits=16)
if(any(abs(good- z[,1])* c(1,1,1,1e-7,1,1e-6) > 1e-7))
      stop("test expo3 sane2.expo3 FAILED")


z <- apply(dfsane1.expo3, 2, summary)
print(z)
print(z[,1], digits=16)
good <- c(
    0,          0,        0,    1.767e-01, 4.450e-07, 2.236e+00)  
#3.357e-08 8.576e-08  9.981e-08 1.767e-01  4.450e-07  2.236e+00 # Linux 32
#3.389e-08 8.126e-08  9.645e-08 1.767e-01  3.110e-07  2.236e+00 # Linux 64
print(good - z[,1], digits=16)
if(any(abs(good  - z[,1])* c(1,1,1,1e-5,1e-2,1e-5) >  1e-7))
     stop("test expo3 dfsane1.expo3 FAILED")


z <- apply(dfsane2.expo3, 2, summary)
print(z)
print(z[,1], digits=16)
good <- c(
    0,          0,        0,      0,          0,        0)
#3.433e-08, 7.543e-08, 8.607e-08, 8.431e-08, 9.332e-08, 2.111e-07 ) # Linux 32
#3.408e-08, 7.616e-08, 8.545e-08, 8.240e-08, 9.247e-08, 9.996e-08 ) # Linux 64
#3.242e-08, 7.543e-08, 8.449e-08, 8.175e-08, 9.261e-08, 1.396e-08   # CRAN Win-builder
print(abs(good  - z[,1])* c(1,1,1,1,1,1e-1), digits=16)
if(any(abs(good - z[,1])* c(1,1,1,1,1,1e-1) > 1e-7))
     stop("test expo3 dfsane2.expo3 FAILED")


#expo3.results <- list(dfsane1=dfsane1.expo3, dfsane2=dfsane2.expo3, sane1=sane1.expo3, sane2=sane2.expo3) 
#dput(expo3.results, file="h:/bb/package/expo3.results")

est1 <- density(log10(sane1.expo3[,1]))
est2 <- density(log10(dfsane1.expo3[,1]))
est3 <- density(log10(sane2.expo3[,1]))
est4 <- density(log10(dfsane2.expo3[,1]))
ymin <- min(est1$y,est2$y,est3$y,est4$y)
ymax <- max(est1$y,est2$y,est3$y,est4$y)
plot(est1,type="l",xlab="log10 (residual)", ylab="kernel density estimate",ylim=c(ymin,ymax),main=" ")
lines(est2,col=2)
lines(est3, lty=2)
lines(est4,col=2, lty=2)
#abline(v=-7, col=4, lty=2, lwd=2)
#legend(locator(1), legend=c("SANE-1","DFSANE-1","SANE-2","DFSANE-2"), lty=c(1,1,2,2),col=c(1,2,1,2))
