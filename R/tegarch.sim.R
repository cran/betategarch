tegarch.sim <-
function(n, omega=0, phi1=0.95, kappa1=0.05,
  kappastar=0, df=10, skew=1, lambda.initial=NULL,
  verbose=FALSE)
{
lambda <- rep(0,n) #lambda
lambdadagg <- rep(0,n) #lambda dagger
if(!is.null(lambda.initial)) lambdadagg[1] <- lambda.initial[2]

epsilon <- rst(n, df=df, skew=skew)
epsilon2 <- epsilon^2
mueps <- st.mean(df=df, skew=skew)
eps2muepseps <- epsilon2 - mueps*epsilon
signeps <- sign(epsilon)
u <- (df+1)*eps2muepseps/(df*skew^(2*signeps) + epsilon2) - 1
signnegyupluss1 <- sign(mueps-epsilon)*(u+1)

#recursion:
fn <- function(i){
  lambdadagg[i+1] <<- phi1*lambdadagg[i] + kappa1*u[i] + kappastar*signnegyupluss1[i]
}
indx <- 1:I(n-1)
lambda1long <- sapply(indx,fn)
lambda <- omega + lambdadagg

#output:
if(verbose){
  u[n] <- NA
  sigma <- exp(lambda)
  epsilon <- epsilon - mueps
  y <- sigma*epsilon
  result <- cbind(y, sigma, lambda, lambdadagg, u, epsilon)
}else{
  sigma <- exp(lambda)
  result <- sigma*(epsilon-mueps)
}
return(result)
}
