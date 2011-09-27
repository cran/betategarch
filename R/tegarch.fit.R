tegarch.fit <-
function(y, par, lambda.initial=NULL, c.code=FALSE){
  tegarch.recursion(y, delta=par[1], phi1=par[2], kappa1=par[3],
  kappa1star=par[4], df=par[5], verbose=TRUE, lambda.initial=lambda.initial,
  c.code=c.code)
}

