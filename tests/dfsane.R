options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


##########
cat("BB test dfsane ...\n")

expo1 <- function(x) {
    	n <- length(x)
    	f <- rep(NA, n)
    	f[1] <- exp(x[1] - 1) - 1
    	f[2:n] <- (2:n) * (exp(x[2:n] - 1) - x[2:n])
    	f
    	}

ans <- dfsane(par=runif(100), fn=expo1)
 
z <- sum(ans$par)
#good   <-  99.9944496188613 #pre 2008-11.1
good   <-   99.9948355416    #2008-11.1
#on Windows 99.9944496188436
#on Linux64 99.9944496188613 pre 2008-11.1
#on Linux32 99.9944496188436 pre 2008-11.1
#on Linux32 99.9948355418059 2008-11.1
#on Linux64 99.9948355414489 2008-11.1
print(z, digits=16)
if(any(abs(good - z) > 1e-9)) stop("BB test dfsane FAILED")
