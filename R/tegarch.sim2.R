tegarch.sim2 <-
function(n, omega=0, phi1=0.95, phi2=0.9,
  kappa1=0.001, kappa2=0.01, kappastar=0.005, df=10, skew=0.98,
  lambda.initial=NULL, verbose=FALSE)
{
lambda <- rep(0,n) #lambda
lambda1dagg <- rep(0,n) #lambda1 dagger
lambda2dagg <- lambda1dagg #lambda2 dagger
if(!is.null(lambda.initial)){
  lambda1dagg[1] <- lambda.initial[2]
  lambda2dagg[1] <- lambda.initial[3]
}

epsilon <- rst(n, df=df, skew=skew)
epsilon2 <- epsilon^2
mueps <- st.mean(df=df, skew=skew)
eps2muepseps <- epsilon2 - mueps*epsilon
signeps <- sign(epsilon)
u <- (df+1)*eps2muepseps/(df*skew^(2*signeps) + epsilon2) - 1
signnegyupluss1 <- sign(mueps-epsilon)*(u+1)

#recursion:
fn <- function(i){
  lambda1dagg[i+1] <<- phi1*lambda1dagg[i] + kappa1*u[i]
  lambda2dagg[i+1] <<- phi2*lambda2dagg[i] + kappa2*u[i] + kappastar*signnegyupluss1[i]
}
indx <- 1:I(n-1); lambda1long <- sapply(indx,fn)
lambda <- omega + lambda1dagg + lambda2dagg
if(!is.null(lambda.initial)) lambda[1] <- lambda.initial[1]

#output:
if(verbose){
  u[n] <- NA
  sigma <- exp(lambda)
  epsilon <- epsilon - mueps
  y <- sigma*epsilon
  result <- cbind(y, sigma, lambda, lambda1dagg, lambda2dagg, u, epsilon)
}else{
  sigma <- exp(lambda)
  epsilon <- epsilon - mueps
  result <- sigma*epsilon
}
return(result)
}
