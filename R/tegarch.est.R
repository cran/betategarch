tegarch.est <-
function(y, initial.values=c(0.001,0.9,0.02,0.01,10),
  lower=c(-Inf,-0.999999999,-Inf,-Inf,2.000001),
  upper=c(+Inf,0.999999999,Inf,Inf,Inf), lambda.initial=NULL,
  compute.hessian=TRUE, c.code=FALSE, na.replace=NA, verbose=TRUE, ...)
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
    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1,
      length(initial.values)), ndeps = rep.int(0.001, length(initial.values)),
      maxit = 100L, abstol = -Inf, reltol = sqrt(.Machine$double.eps),
      alpha = 1, beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
      factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
    hess <- .Internal(optimhess(est$par, fn1, NULL, con))
    est$hessian.numerical <- 0.5 * (hess + t(hess))
    est$hessian.numerical <- -est$hessian.numerical
  }

if(verbose){
  names(est$par) <- c("delta", "phi1", "kappa1", "kappa1star", "df")
}

return(est)
}

