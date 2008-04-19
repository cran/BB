
options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("BB test broydt.f ...\n")

broydt.f <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
sum(f*f)
}

p0 <- rnorm(100, sd=1)
system.time(ans.spg <- spg(par=p0, fn=broydt.f))[1]
system.time(ans.opt <- optim(par=p0, fn=broydt.f, method="L-BFGS-B"))[1]
 
z <- sum(ans.spg$par)
good   <-   98.55031219348329
#on Windows 98.513096595265
#on Linux64 98.55031219348329
#on Linux32 98.6221556927617
#on CRAN    98.6625515502629
print(z, digits=16)
if(any(abs(good - z) > 5e-1)) stop("BB test broydt.f a FAILED")
 
z <- sum(ans.opt$par)
good   <-   111.5078705487698
#on Windows 111.5442515847844
#on Linux64 111.5078705487698
#on Linux32 111.5128927193081
print(z, digits=16)
if(any(abs(good - z) > 5e-1)) stop("BB test broydt.f b FAILED")
