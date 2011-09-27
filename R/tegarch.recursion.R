tegarch.recursion <-
function(y, delta=0.01, phi1=0.95, kappa1=0.05,
  kappa1star=0, df=10, lambda.initial=NULL, c.code=FALSE, verbose=FALSE)
{
y2 <- y^2
kappa1starsignnegy <- kappa1star*sign(-y)
dfpluss1y2 <- (df+1)*y2
iN <- length(y)
u <- rep(0,iN)
lambda1 <- rep(0,iN)
if(is.null(lambda.initial)){
  lambda1[1] <- if(abs(phi1)==1){ 0 }else{ delta/(1-phi1) }
}else{ lambda1[1] <- lambda.initial }

if(c.code){
  tmp <- .tegarch.recursion(as.integer(iN), as.numeric(delta),
    as.numeric(phi1), as.numeric(kappa1), as.numeric(df), as.numeric(y2),
    as.numeric(dfpluss1y2), as.numeric(kappa1starsignnegy),
    as.numeric(lambda1), as.numeric(u))
  u <- tmp$u
  u[iN] <- NA
  lambda1 <- tmp$lambda1
}else{
  fn <- function(i){
    u[i] <<- dfpluss1y2[i]/(df*exp(2*lambda1[i]) + y2[i]) - 1
    lambda1[i+1] <<- delta + phi1*lambda1[i] + kappa1*u[i] + kappa1starsignnegy[i]*(u[i]+1)
  }
  indx <- 1:I(iN-1)
  tmp <- sapply(indx,fn)
}

#output:
if(verbose){
  sigma <- exp(lambda1)
  epsilon <- y/sigma
  result <- cbind(y,sigma,lambda1,u,epsilon)
}else{ result <- lambda1 }

return(result)
} #end tegarch.recursion

