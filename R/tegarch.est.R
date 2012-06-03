tegarch.est <-
function(y, initial.values=c(0.001,0.9,0.02,0.01,10),
  lower=c(-Inf,-0.999999999,-Inf,-Inf,2.000001),
  upper=c(+Inf,0.999999999,Inf,Inf,Inf), lambda.initial=NULL,
  compute.hessian=TRUE, c.code=TRUE, na.replace=NA, verbose=TRUE, ...)
{
  objective.f <- function(p, x=y, logl.scale=-1){
    f <- logl.scale*tegarch.logl(x, p[1],p[2],p[3],p[4],p[5], lower=lower,
    upper=upper, lambda.initial=lambda.initial, na.replace=na.replace,
    c.code=c.code); f}
  est <- nlminb(initial.values, objective.f, lower=lower, upper=upper,
    ..., x=y)
  est$objective <- -est$objective

  if(compute.hessian){
    fn1 <- function(initial.values) objective.f(initial.values, logl.scale=-1)
    est$hessian.numerical <- -optimHess(est$par, fn1)
  }

if(verbose){
  names(est$par) <- c("delta", "phi1", "kappa1", "kappa1star", "df")
}

return(est)
}
