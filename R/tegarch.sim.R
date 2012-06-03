tegarch.sim <-
function(n, delta=0, phi1=0.95, kappa1=0.05,
  kappa1star=0, df=10, lambda.initial=NULL, verbose=FALSE)
{
  lambda1 <- rep(NA,n)
  if(is.null(lambda.initial)){
    lambda1[1] <- if(abs(phi1)==1){ 0 }else{ delta/(1-phi1) }
  }else{ lambda1[1] <- lambda.initial }

  epsilon <- rt(n, df=df)
  epsilon2 <- epsilon^2
  signy <- sign(epsilon)
  u <- (df+1)*epsilon2/(df + epsilon2) - 1
  uterm <- delta + kappa1*u
  if(kappa1star != 0){
    uterm <- uterm + kappa1star*I(-1)*signy*(u+1)
  }

  #recursion:
  fn <- function(i){lambda1[i+1] <<- phi1*lambda1[i] + uterm[i]}
  indx <- 1:I(n-1); lambda1long <- sapply(indx,fn)

  #output:
  if(verbose){
    u[n] <- NA
    sigma <- exp(lambda1)
    y <- sigma*epsilon
    result <- cbind(y, sigma, lambda1, u, epsilon)
  }else{
    sigma <- exp(lambda1)
    result <- sigma*epsilon
  }
  return(result)
}
