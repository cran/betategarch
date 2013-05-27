tegarch.est <-
function(y, asym=TRUE, skew=TRUE,
  components=1, initial.values=NULL, lower=NULL, upper=NULL,
  compute.hessian=FALSE, lambda.initial=NULL, c.code=TRUE,
  logl.penalty=-1e+100, verbose=TRUE, aux=NULL, ...)
{
#if(verbose){
#  y <- as.zoo(y)
#  y <- na.trim(y)
#  y.index <- index(y)
#  y <- coredata(y)
#}
aux$asym <- asym
aux$skew <- skew
aux$iN <- length(y)
aux$signnegy <- sign(-y)
aux$u <- rep(0,aux$iN)

if(components==1){
  if(is.null(initial.values)){
    initial.values <- c(0.02,0.95,0.05,0.01,10,0.98)
  }
  if(is.null(lower)){
    lower <- c(-Inf,-1+.Machine$double.eps,-Inf,-Inf,
      2+.Machine$double.eps,.Machine$double.eps)
  }
  if(is.null(upper)){
    upper <- c(Inf,1-.Machine$double.eps,Inf,Inf,Inf,Inf)
  }
  if(!aux$skew){
    initial.values <- initial.values[-6]
    lower <- lower[-6]
    upper <- upper[-6]
  }
  if(!aux$asym){
    initial.values <- initial.values[-4]
    lower <- lower[-4]
    upper <- upper[-4]
  }
  if(is.null(logl.penalty)){
    logl.penalty <- tegarch.logl(y, initial.values,
      lower=lower, upper=upper, lambda.initial=lambda.initial,
      logl.penalty=-1e+100, c.code=c.code, aux=aux)
  }
  objective.f <- function(pars, x=y){f <- -tegarch.logl(x,
    pars, lower=lower, upper=upper,
    lambda.initial=lambda.initial, logl.penalty=logl.penalty,
    c.code=c.code, aux=aux); f}
}else{
  if(is.null(initial.values)){
    initial.values <- c(0.02,0.95,0.9,0.001,0.01,0.005,10,0.98)
  }
  if(is.null(lower)){
    lower <- c(-Inf,-1+.Machine$double.eps,
      -1+.Machine$double.eps,-Inf,-Inf,-Inf,
      2+.Machine$double.eps,.Machine$double.eps)
  }
  if(is.null(upper)){
  upper <- c(Inf,1-.Machine$double.eps,
    1-.Machine$double.eps,Inf,Inf,Inf,Inf,Inf)
  }
  if(!aux$skew){
    initial.values <- initial.values[-8]
    lower <- lower[-8]
    upper <- upper[-8]
  }
  if(!aux$asym){
    asym=TRUE
  }
  if(is.null(logl.penalty)){
    logl.penalty <- tegarch.logl2(y, initial.values,
      lower=lower, upper=upper, lambda.initial=lambda.initial,
      logl.penalty=-1e+100, c.code=c.code, aux=aux)
  }
  objective.f <- function(pars, x=y){f <- -tegarch.logl2(x,
    pars, lower=lower, upper=upper, lambda.initial=lambda.initial,
    logl.penalty=logl.penalty, c.code=c.code, aux=aux); f}
}

#estimate:
#est <- nlminb(initial.values, objective.f, lower=lower,
#  upper=upper, x=y)
est <- nlminb(initial.values, objective.f, lower=lower,
  upper=upper, x=y, ...)
est$objective <- -est$objective

#compute Hessian:
if(compute.hessian){
  hessian.numerical <- -optimHess(est$par, objective.f)
  est <- c(list(hessian.numerical=hessian.numerical), est)
}

#type of model:
model <- c(components, asym, skew)
names(model) <- c("components", "asym", "skew")
est <- c(list(model=model), est)

if(verbose){
  if(components==1){
    parnames <- c("omega", "phi1", "kappa1", "kappastar", "df", "skew")
    if(!aux$skew){ parnames <- parnames[-6] }
    if(!aux$asym){ parnames <- parnames[-4] }
  }else{
    parnames <- c("omega", "phi1", "phi2", "kappa1", "kappa2",
      "kappastar", "df", "skew")
    if(!aux$skew){ parnames <- parnames[-8] }
    if(!aux$asym){
      est$NOTE <- "2 comp spec without leverage not available"
    }
  }
  names(est$par) <- parnames
  if(compute.hessian){
    colnames(est$hessian.numerical) <- parnames
    rownames(est$hessian.numerical) <- parnames
  }
  names(initial.values) <- parnames
  names(lower) <- parnames
  names(upper) <- parnames
  est <- c(list(date=date(), initial.values=initial.values,
    lower=lower, upper=upper), est)
}

return(est)
}
