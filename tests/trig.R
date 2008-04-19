options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


##########
cat("BB test trig.f ...\n")

trig.f <- function(x){
n <- length(x)
i <- 1:n
f <- n - sum(cos(x)) + i*(1 - cos(x)) - sin(x) 
sum(f*f)
}

p0 <- rnorm(500,sd=5)
system.time(ans.spg <- spg(par=p0, fn=trig.f, control=list(maxit=2500)))[1]
system.time(ans.opt <- optim(par=p0, fn=trig.f, method="L-BFGS-B", control=list(lmm=10)))[1]
 
z <- sum(ans.spg$par)
good   <-   126.5285048574777
#on Windows 126.5285707910017
#on Windows 126.5285707910017 on Uwe's win-builder
#on Linux64 126.5285048574777
#on Linux32 126.5285013419033
print(z, digits=16)
if(any(abs(good - z) > 1e-4)) stop("BB test trig.f a FAILED")
 
z <- sum(ans.opt$par)
good   <-   126.6745362827593
#on Windows 126.6745362708396 
#on Windows 126.670123957483  on Uwe's win-builder
#on Linux64 126.6745362827593
#on Linux32 126.6745362418508
print(z, digits=16)
if(any(abs(good - z) > 1e-2)) stop("BB test trig.f b FAILED")
