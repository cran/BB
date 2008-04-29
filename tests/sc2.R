options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

#########################################################################################
cat("BB test sc2 ...\n")

sc2.f <- function(x){
n <- length(x)
vec <- 1:n
sum(vec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
n <- length(x)
vec <- 1:n
vec * (exp(x) - 1) / 10
}

neg.sc2.f <- function(x){
n <- length(x)
vec <- 1:n
-sum(vec * (exp(x) - x)) / 10
}

neg.sc2.g <- function(x){
n <- length(x)
vec <- 1:n
-vec * (exp(x) - 1) / 10
}

p0 <- runif(500,min=-1, max=1)
system.time(ans.spg <- spg(par=p0, fn=sc2.f, control=list(maxit=2500)))[1]

z <- sum(ans.spg$par)
good   <-    2.862543216705235e-05
#on Windows -0.0002158022390025393
print(z, digits=16)
if(any(abs(good - z) > 1e-3)) stop("BB test sc2 a FAILED")

system.time(neg.ans.spg <- spg(par=p0, fn=neg.sc2.f, 
              control=list(maxit=2500, maximize=TRUE)))[1]

z <- sum(neg.ans.spg$par)
good   <-    2.862543216705235e-05
print(z, digits=16)
if(any(abs(good - z) > 1e-3)) stop("BB test neg sc2 a FAILED")

system.time(ans.spg <- spg(par=p0, fn=sc2.f, gr=sc2.g,
   control=list(maxit=2500)))[1]

z <- sum(ans.spg$par)
good <- 2.565413040899874e-06
#on Linux64 2.565413040899874e-06 (mfacl2)
#on Linux64 6.677493403589264e-05
print(z, digits=16)
if(any(abs(good - z) >  1e-4)) stop("BB test sc2 b FAILED")

system.time(neg.ans.spg <- spg(par=p0, fn=neg.sc2.f, gr=neg.sc2.g,
   control=list(maxit=2500, maximize=TRUE)))[1]

z <- sum(neg.ans.spg$par)
good <- 2.565413040899874e-06
print(z, digits=16)
if(any(abs(good - z) >  1e-4)) stop("BB test neg.sc2 b FAILED")

system.time(ans.opt <- optim(par=p0, fn=sc2.f, method="L-BFGS-B"))[1]

z <- sum(ans.opt$par)
good   <-   0.02209066162550582
#on Windows 0.02209186415471651
#on Linux64 0.02209066162550582
#on Linux32 0.0220908989551237
print(z, digits=16)
if(any(abs(good - z) > 1e-4)) stop("BB test sc2 c FAILED")

system.time(ans.opt <- optim(par=p0, fn=sc2.f, gr=sc2.g, method="L-BFGS-B"))[1]

z <- sum(ans.opt$par)
good   <-   0.02200130759852783
#on Windows 0.02200130758634488
#on Linux64 0.02200130759852783
#on Linux32 0.02200130759779934
print(z, digits=16)
if(any(abs(good - z) > 1e-9)) stop("BB test sc2 d FAILED")
