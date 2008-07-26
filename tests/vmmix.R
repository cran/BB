options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

###############################################################
cat("BB test vmmix ...\n")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

dvm <- function (theta, mu, kappa) 
{
    1/(2 * pi * besselI(x = kappa, nu = 0, expon.scaled = TRUE)) * 
        (exp(cos(theta - mu) - 1))^kappa
}

##
rmixedvm <- function (n, mu1, mu2, kappa1, kappa2, p) {
temp <- runif(n)
n1 <- sum(temp <= p)
y <- c(rvm(n1,mu1,kappa1),rvm(n-n1,mu2,kappa2))
return(y)
}

##
rvm <- function (n, mean, k) 
{
    vm <- c(1:n)
    a <- 1 + (1 + 4 * (k^2))^0.5
    b <- (a - (2 * a)^0.5)/(2 * k)
    r <- (1 + b^2)/(2 * b)
    obs <- 1
    while (obs <= n) {
        U1 <- runif(1, 0, 1)
        z <- cos(pi * U1)
        f <- (1 + r * z)/(r + z)
        c <- k * (r - f)
        U2 <- runif(1, 0, 1)
        if (c * (2 - c) - U2 > 0) {
            U3 <- runif(1, 0, 1)
            vm[obs] <- sign(U3 - 0.5) * acos(f) + mean
            vm[obs] <- vm[obs]%%(2 * pi)
            obs <- obs + 1
        }
        else {
            if (log(c/U2) + 1 - c >= 0) {
                U3 <- runif(1, 0, 1)
                vm[obs] <- sign(U3 - 0.5) * acos(f) + mean
                vm[obs] <- vm[obs]%%(2 * pi)
                obs <- obs + 1
            }
        }
    }
    vm
}

#
vmmix.loglik <- function(x,y){
p <- x
- sum(log(p[5]*dvm(y,p[1],p[2])+(1-p[5])*dvm(y,p[3],p[4])))
}

y <- rmixedvm(n=500, mu1=pi/2, mu2=3*pi/2, kappa1=1.9, kappa2=2.2, p=0.67)
p <- c(pi/4,2,pi,1,0.5)

lo <- rep(0.001,5)
hi <- c(Inf, Inf, Inf, Inf, 0.999)

p <- c(runif(5,c(0,0.1,0,0.1,0.2),c(2*pi,5,2*pi,5,0.8)))
system.time(ans.spg <- spg(par=p, fn=vmmix.loglik, y=y, lower=lo, upper=hi, 
   control=list(maxit=2500, gtol=1.e-06, M=20, trace=T)))[1]
ans.opt <- optim(par=p, fn=vmmix.loglik, y=y, method="L-BFGS-B", lower=lo, upper=hi)

# this should give the same reult. It just tests passing project arg
Userproject <- function(x, lower, upper, ...) {
       x[x < lower] <- lower[x < lower]
       x[x > upper] <- upper[x > upper]
       return(x)
       }

ans.spg2 <- spg(par=p, fn=vmmix.loglik, y=y, lower=lo, upper=hi,
   project=Userproject,
   control=list(maxit=2500, gtol=1.e-06, M=20, trace=T))[1]


z <- sum(ans.spg$par)
good   <-   11.28100000000000
#on Windows 11.28096954777041
#on Linux64 11.28094286243805
#on Linux32 11.28105871638063
#on Linux32 11.28061301609979
print(z, digits=16)
if(any(abs(good - z) > 1e-3)) stop("BB test vmmix.loglik a1 FAILED")


z <- sum(ans.spg2$par)
print(z, digits=16)
if(any(abs(good - z) > 1e-3)) stop("BB test vmmix.loglik a2 FAILED")
 
z <- sum(ans.opt$par)
good   <-   11.28032840876373
#on Windows 11.28032868302692
#on Linux64 11.28032840876373
#on Linux32 11.28032868197765
print(z, digits=16)
if(any(abs(good - z) > 1e-4)) stop("BB test vmmix.loglik  b FAILED")

# previously had this, but it requires numDeriv
#gs <- grad(ans.spg$par, func=vmmix.loglik, y=y)
#go <- grad(ans.opt$par, func=vmmix.loglik, y=y)
#
#z <- sum(gs)
#good	<-   -0.0001803984142419232
##on Windows -0.00016745341957699
##on Linux64 -0.0001803984142419232
##on Linux32 -0.000008890766665183735
#print(z, digits=16)
#if(any(abs(good - z) > 1e-3)) stop("BB test vmmix.loglik a FAILED")
#
#z <- sum(go )
#good	<-   -0.02992184037618598
##on Windows -0.02989100666129558
##on Linux64 -0.02992184037618598
##on Linux32 -0.029891113567616
#print(z, digits=16)
#if(any(abs(good - z) > 1e-3)) stop("BB test vmmix.loglik  b FAILED")
