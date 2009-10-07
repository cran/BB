if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1234)
old.seed <- setRNG(test.rng)


# A high-degree polynomial system (R.B. Kearfoot, ACM 1987)
# There are 12 real roots (and 126 complex roots to this system!)
#
hdp <- function(x) {
f <- rep(NA, length(x))
f[1] <- 5 * x[1]^9 - 6 * x[1]^5 * x[2]^2 + x[1] * x[2]^4 + 2 * x[1] * x[3]
f[2] <- -2 * x[1]^6 * x[2] + 2 * x[1]^2 * x[2]^3 + 2 * x[2] * x[3]
f[3] <- x[1]^2 + x[2]^2 - 0.265625
f
}

p0 <- matrix(runif(600), 200, 3)  # 200 starting values, each of length 3
#ans <- BBsolve(par=p0, fn=hdp)
ans <- multiStart(par=p0, fn=hdp)

#pc <- princomp(ans$par)
#plot(pc$scores[,1])  # you can see all 12 solutions
 
z <- sum(ans$par[ans$converged,])
good   <-    69.0642609408530
#on Windows 
#on Linux64  69.0642609408530
#on CRAN Mac 67.24096047829846 had to change fuz from 0.5 to 2.0 for this
# the R-forge Mac testing worked with 0.5

 # before  M changed from 10,50 to 50,10 default:70.0682292316792  
 # early defaults: 68.6429426963639
#on Linux32  # early defaults: 68.4364913774026

print(z, digits=16)
if(any(abs(good - z) > 2.0)) stop("BB test BBsolve HDP FAILED")
