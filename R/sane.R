 sane <- function(par, fn, method=1, control=list(), ...) {

  # control defaults
  ctrl <- list(M=10, maxit=1500, tol=1e-07, trace=TRUE, triter=10) 
  namc <- names(control)
  if (! all(namc %in% names(ctrl)) )
     stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])     

  ctrl[namc ] <- control
  M	 <- ctrl$M
  maxit  <- ctrl$maxit
  tol	 <- ctrl$tol
  trace  <- ctrl$trace
  triter <- ctrl$triter

  fargs <- list(...)
  #######################  local function
  #  non-monotone line search of Grippo
  lineSearch <- function(x, fn, F, fval, dg, M, lastfv, sgn, lambda, 
                     fcnt, bl, fargs) {
  maxbl  <- 100
  gamma  <- 1.e-04
  sigma1 <- 0.1
  sigma2 <- 0.5
  cbl    <- 0
  fmax   <- max(lastfv)
  gpd    <- -2 * abs(dg)
  
  #    line search Main loop
  while (cbl < maxbl) {

    xnew <-  x + lambda* sgn * F

    Fnew <- try(do.call("fn", append(list(xnew) , fargs )))
    fcnt = fcnt + 1

    if (class(Fnew) == "try-error" | any(is.nan(Fnew)) )
         return(list(xnew=NA, Fnew=NA, fcnt=fcnt, bl=bl, lsflag=1, fune=NA))
      else fune <- sum(Fnew * Fnew)

    if (fune <= (fmax + lambda*gpd*gamma)) {
      if (cbl >= 1) bl <- bl + 1
           return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, lambda=lambda, bl=bl, 
	               lsflag=0, fune=fune))
      } else {
        #     Quadratic interpolation
        lamc <- -(gpd*lambda^2) / (2 * (fune - fval - lambda*gpd))
        c1 <- sigma1 * lambda
        c2 <- sigma2 * lambda
        if (lamc < c1) lambda <- c1
        else if (lamc > c2)lambda <- c2
        else lambda <- lamc
        cbl <- cbl + 1
        }
    }
  return(list(xnew=NA, Fnew=NA, fcnt=fcnt, lambda=NA, bl=bl, lsflag=2, fune=NA))
  }

  ##########################################

  #     Initialization
  n      <- length(par)
  fcnt   <- iter <- bl <- 0
  alfa   <- 1
  eps    <- 1.0e-08
  h      <- 1.0e-07
  lastfv <- rep(0, M)
  stagn  <- FALSE

  F <- try (fn(par, ...))

  if (class(F) == "try-error" )
     stop(" Failure in initial functional evaluation.", F) 
  else if (any(is.nan(F), is.infinite(F), is.na(F)) )
     stop(" Failure in initial functional evaluation. Try another initial value.") 
  
  F0 <- normF <- sqrt(sum(F * F))

  if (trace) cat("Iteration: ", 0, " ||F(x0)||: ", F0, "\n")

  pbest <- par
  normF.best <- normF
  lastfv[1] <- normF^2
  	  
  ######################################
  #     MAIN LOOP:  Iteration begins
  ######################################
  while (normF/sqrt(n) > tol & iter <= maxit & !stagn) {

    # Calculate the gradient of the merit function ||F(X)||

    Fa <- try (fn(par + h*F, ...)) 
    fcnt <- fcnt + 1
        
    if (class(Fa) == "try-error" | any(is.nan(Fa)) ) {
        flag <- 1
        break
        } 

    dg <- (sum(F * Fa) - normF^2) / h
    if (abs(dg/normF^2) < eps | is.nan(dg) | is.infinite(dg) ) {
        flag <- 3
        break
        }

    # Control of steplength
    if ((alfa <= eps) | (alfa >= 1/eps)) 
        alfa <- if (normF > 1) 1 
                else if (normF >= 1e-05 & normF <= 1)   normF
                else if (normF < 1e-05)  1e-05
 
    sgn <- if (dg > 0) -1 else  1

    lambda <- 1/alfa

    #  non-monotone line search of Grippo
    ls.ret <-  lineSearch(x=par, fn=fn, F=F, fval=normF^2, dg=dg, M=M, 
                lastfv=lastfv, sgn, lambda, fcnt, bl, fargs)

    fcnt <- ls.ret$fcnt
    bl <- ls.ret$bl
    flag <- ls.ret$lsflag

    if (flag > 0) break 

    lambda <- ls.ret$lambda
    Fnew   <- ls.ret$Fnew
    pnew   <- ls.ret$xnew
    fune   <- ls.ret$fune

    #    Calculate new steplength: alfa

    alfa <- 
       if (method==1) sum(F*(F - Fnew)) / (lambda * sum(F*F)) else
       if (method==2) sum((F - Fnew)^2) / (lambda * sum(F*(F - Fnew))) else
       if (method==3) sign(sum(F*(F - Fnew))) * sqrt(sum((F - Fnew)^2) / (lambda^2 * sum(F*F)))

    if (is.nan(alfa)) alfa <- 1.e-08
    
    par <- pnew
    F   <- Fnew
    fun <- fune
    normF <- sqrt(fun)

    if (normF < normF.best) {
      pbest <- par
      normF.best <- normF
      }

    iter <- iter + 1
    lastfv[ 1 + iter %% M] <- fun
    
    if (trace && (iter%%triter == 0))
       cat("\n iteration: ",iter, " ||F(xn)|| =  ", normF, "\n")

    }  # End of main loop

  if (flag==0) {
     if (normF/sqrt(n) <= tol) 
                  conv <- list(type=0, message="Successful convergence") 
     if (iter > maxit) 
                  conv <- list(type=1, message="Maximum number of iterations exceeded")
     } else {
     par <- pbest
     normF <- normF.best
     if (flag==1) conv <- list(type=2, message="Failure: Error in function evaluation")
     if (flag==2) conv <- list(type=3, message="Failure: Maximum limit on steplength reductions exceeded")
     if (flag==3) conv <- list(type=4, message="Failure: Anomalous iteration")
     }
  
  return(list(par = par, 
     residual = normF / sqrt(length(par)), 
     fn.reduction = F0 - normF, 
     feval=fcnt, iter=iter, convergence=conv$type, message=conv$message)) 
  }


