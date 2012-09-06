.tegarch.recursion <- function (iN, delta, phi1, kappa1, df, y2, dfpluss1y2, kappa1starsignnegy, 
    lambda1, u) 
.C("file1de03a27a74", iN = as.integer(iN), delta = as.double(delta), 
    phi1 = as.double(phi1), kappa1 = as.double(kappa1), df = as.double(df), 
    y2 = as.double(y2), dfpluss1y2 = as.double(dfpluss1y2), kappa1starsignnegy = as.double(kappa1starsignnegy), 
    lambda1 = as.double(lambda1), u = as.double(u), PACKAGE = "betategarch")