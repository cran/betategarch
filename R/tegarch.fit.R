tegarch.fit <-
function(y, result, lambda.initial=NULL)
{
if(result$model["components"]==1){
  if(result$model["skew"]==0){
    pars <- c(result$par, 1)
  }else{ pars <- result$par }
  if(length(pars) < 6){
    pars <- c(pars[1:3], 0, pars[4:5])
  }
  out <- tegarch.recursion(y, omega=pars[1], phi1=pars[2],
    kappa1=pars[3], kappastar=pars[4], df=pars[5],
    skew=pars[6], lambda.initial=lambda.initial,
    c.code=TRUE, verbose=TRUE, aux=NULL)
}else{
  if(result$model["skew"]==0){
    pars <- c(result$par, 1)
  }else{ pars <- result$par }
  out <- tegarch.recursion2(y, omega=pars[1], phi1=pars[2],
    phi2=pars[3], kappa1=pars[4], kappa2=pars[5],
    kappastar=pars[6], df=pars[7], skew=pars[8],
    lambda.initial=lambda.initial, c.code=TRUE,
    verbose=TRUE, aux=NULL)
}

return(out)
}
