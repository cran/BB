# test results in these files are indicated for
# Windows WP Professional v5.1 SP2
# Linux32 Ubuntu 7.10 desktop 32 bit kernel 2.6.22-14-generic on Intel Pentium
# Linux64 Gentoo   64 bit kernel 2.6.17-gentoo-r8 on AMD 64 X2

options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("BB test brown.f ...\n")

brown.f <- function(x) {
p <- x
n <- length(p)
odd <- seq(1,n,by=2)
even <- seq(2,n,by=2)
sum((p[odd]^2)^(p[even]^2 + 1) + (p[even]^2)^(p[odd]^2 + 1))
}

p0 <- rnorm(500,sd=2) # this set fails in optim, so
p0 <- rnorm(500,sd=2)
system.time(ans.spg <- spg(par=p0, fn=brown.f, control=list(maxit=2500)))[1]

z <- sum(ans.spg $par)
good  <-    -3.336592920523865e-05 
#on Windows -3.337359033283664e-05
#on Linux64 -3.383384527412175e-05
#on Linux32 -3.336592920523865e-05
print(z, digits=16)
if(any(abs(good - z) > 1e-6)) stop("BB test brown.f a FAILED")

system.time(ans.opt <- optim(par=p0, fn=brown.f, method="L-BFGS-B"))[1]

z <- sum(ans.opt $par)
good   <-   0.00898911426823444
#on Windows 0.00898911689533186
#on Linux64 0.00898911426823444
#on Linux32 0.00898911750186776
print(z, digits=16)
if(any(abs(good - z) > 1e-7)) stop("BB test brown.f b FAILED")
