tegarch.logl <-
function(y, delta=0, phi1=0.8, kappa1=0.1,
  kappa1star=0.01, df=10, lower=c(-Inf,-0.999999999,-Inf,-Inf,2.1),
  upper=c(Inf,0.999999999,Inf,Inf,Inf), lambda.initial=NULL, c.code=TRUE,
  na.replace=rep(NA,5))
{
if(is.na(delta)) delta <- na.replace[1]
if(delta <= lower[1]) delta <- lower[1]
if(delta >= upper[1]) delta <- upper[1]
if(is.na(phi1)) phi1 <- na.replace[2]
if(phi1 <= lower[2]) phi1 <- lower[2]
if(phi1 >= upper[2]) phi1 <- upper[2]
if(is.na(kappa1)) kappa1 <- na.replace[3]
if(kappa1 <= lower[3]) kappa1 <- lower[3]
if(kappa1 >= upper[3]) kappa1 <- upper[3]
if(is.na(kappa1star)) kappa1star <- na.replace[4]
if(kappa1star <= lower[4]) kappa1star <- lower[4]
if(kappa1star >= upper[4]) kappa1star <- upper[4]
if(is.na(df)) df <- na.replace[5]
if(df <= lower[5]) df <- lower[5]
if(df >= upper[5]) df <- upper[5]

lambda1 <- tegarch.recursion(y, delta=delta, phi1=phi1, kappa1=kappa1,
  kappa1star=kappa1star, df=df, lambda.initial=lambda.initial,
  c.code=c.code)

iN <- length(y)
y2 <- y^2
denom.term <- df*exp(2*lambda1)

const1 <- lgamma((df+1)/2) - lgamma(df/2)  - log(pi*df)/2
term1 <- lambda1
term2 <- (df+1)*log(1 + (y2/denom.term))/2

logl <- iN*const1 - sum(term1) - sum(term2)
if(is.na(logl) || abs(logl) == Inf) logl <- -10e+100
return(logl)
}

