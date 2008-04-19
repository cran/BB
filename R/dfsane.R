dfsane <- function(par, fn, method=2, control=list(), ...) {

  # control defaults
  ctrl <- list(M=10, maxit=1500, tol=1e-07, trace=TRUE, triter=10) 
  namc <- names(control)
  if (! all(namc %in% names(ctrl)) )
     stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])     

  ctrl[namc ] <- control
  M	<- ctrl$M
  maxit <- ctrl$maxit
  tol	<- ctrl$tol
  trace <- ctrl$trace
  triter <- ctrl$triter

  fargs <- list(...)
  ############   local function
  lsm <- function(x, fn, F, fval, alfa, M, lastfv, eta, fcnt, bl, fargs) {
    #  non-monotone line search of Grippo
  
    maxbl <- 100
    gamma <- 1.e-04
    sigma1 <- 0.1
    sigma2 <- 0.5
    lam1 <-  lam2 <- 1
    cbl <- 0
    fmax <- max(lastfv)
    
    #  local function Main loop
    while (cbl < maxbl) {

    	d <- - alfa * F
    	xnew <-  x + lam1 * d

    	Fnew <- try(do.call("fn", append(list(xnew) , fargs )))
    	fcnt = fcnt + 1

    	if (class(Fnew) == "try-error" | any(is.nan(Fnew)) )
    	   return(list(xnew=NA, Fnew=NA, fcnt=fcnt, bl=bl, lsflag=1, fune=NA))
    	fune1 <- sum(Fnew * Fnew)

    	if (fune1 <= (fmax + eta - (lam1^2*gamma*fval))) {
    	    if (cbl >= 1) bl <- bl + 1
    	    return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=0, fune=fune1))
        } 

    xnew <- x - lam2 * d

    Fnew  <- try(do.call("fn", append(list(xnew) , fargs )))
    fcnt = fcnt + 1

    if (class(Fnew) == "try-error" | any(is.nan(Fnew)) )
         return(list(xnew=NA, Fnew=NA, fcnt=fcnt, bl=bl, lsflag=1, fune=NA)) 

    fune2 <- sum(Fnew * Fnew)

    if (fune2 <= (fmax + eta - (lam2^2*gamma*fval))) {
      if (cbl >= 1) bl <- bl + 1
      return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=0, fune=fune2))
      } 
    #     Quadratic interpolation
    lamc <- (2*fval*lam1^2) / (2 * (fune1 + (2*lam1 - 1)*fval))
    c1 <- sigma1 * lam1
    c2 <- sigma2 * lam1
    lam1 <- if (lamc < c1) c1 else if (lamc > c2) c2 else lamc

    lamc <- (2*fval*lam2^2) / (2 * (fune2 + (2*lam2 - 1)*fval))
    c1 <- sigma1 * lam2
    c2 <- sigma2 * lam2
    lam2 <- if (lamc < c1) c1 else if (lamc > c2) c2 else lamc

    cbl <- cbl + 1
    }

    return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=2, fune=fune))
    }

##########################################
#
#     Initialization
    n <- length(par)
    fcnt <- iter <- bl <- kount <- 0
    alfa <- eta <- 1
    ep <- 1.0e-10
    lastfv <- rep(0, M)
      
    F <- try (fn(par, ...))

    if (class(F) == "try-error" ) 
      stop("Failure in initial functional evaluation.", F)
    else if (any(is.nan(F), is.infinite(F), is.na(F)) ) 
      stop("Failure in initial functional evaluation. Try another starting value.")
    
    F0 <- normF <- sqrt(sum(F * F))

    if (trace) cat("Iteration: ", 0, " ||F(x0)||: ", F0, "\n")
    pbest <- par
    normF.best <- normF

    lastfv[1] <- normF^2
            
######################################
#     MAIN LOOP:  Iteration begins
######################################
    while (normF/sqrt(n) > tol & iter <= maxit) {
 
    	# Control of steplength
    	if ((abs(alfa) <= ep) | (abs(alfa) >= 1/ep))
    	    alfa <- if (normF > 1)   1    else 
		    if (normF >= 1e-05 & normF <= 1) 1/normF else 
		    if (normF < 1e-05)   1e05 

    	#  non-monotone line search of Grippo
    	ls.ret <-  lsm(x=par, fn=fn, F=F, fval=normF^2, alfa, M=M,
    		       lastfv=lastfv, eta, fcnt, bl, fargs)
 
    	fcnt <- ls.ret$fcnt
    	bl   <- ls.ret$bl
    	flag <- ls.ret$lsflag

    	if (flag > 0) {
    	  par <- pbest
    	  normF <- normF.best
    	  break
    	  }
 
    	Fnew <- ls.ret$Fnew
    	pnew <- ls.ret$xnew
    	fune <- ls.ret$fune
 
    	#	  Calculate new steplength: alfa
    	#	Three Barzilai-Borwein steplength options.
    	#
    	pF <- sum((pnew - par) * (Fnew - F))
    	pp <- sum((pnew - par)^2)
    	FF <- sum((Fnew - F)^2)

    	alfa <-  (if (method==1)  pp / pF
    	     else if (method==2)  pF / FF
    	     else if (method==3) sign(pF) * sqrt(pp / FF) )

    	if (is.nan(alfa)) alfa <- 1.e-10

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
    	eta <- F0 / (iter + 1)^2
 
    	if (trace && (iter%%triter == 0))
    	    cat("iteration: ",iter, " ||F(xn)|| =  ", normF, "\n")

    	}   # End of main loop

    conv <-
     if (flag==0) {
       if (normF/sqrt(n) <= tol) 
            list(type=0, message="Successful convergence") else 
       if (iter > maxit)
            list(type=1, message="Maximum limit for iterations exceeded")
       }  else   
       if (flag==1)
            list(type=2, message="Failure: Error in function evaluation") else 
       if (flag==2)
            list(type=3, message="Failure: Maximum limit on steplength reductions exceeded")


    return(list(par = par, 
      residual = normF / sqrt(length(par)), 
      fn.reduction = F0 - normF, 
      feval=fcnt, 
      iter=iter, convergence=conv$type, 
      message=conv$message)) }
