#include <R.h>


extern "C" {
  void tegarchrecursion ( int * iN, double * omega, double * phi1, double * kappa1, double * kappastar, double * df, double * skew2, double * dfpluss1, double * mueps, double * y, double * signnegy, double * signarg, double * skewterm, double * lambda, double * lambdadagg, double * u );
}

void tegarchrecursion ( int * iN, double * omega, double * phi1, double * kappa1, double * kappastar, double * df, double * skew2, double * dfpluss1, double * mueps, double * y, double * signnegy, double * signarg, double * skewterm, double * lambda, double * lambdadagg, double * u ) {

  for (int i=0; i < *iN-1; i++) {
    signarg[i] = y[i] + *mueps * exp(lambda[i]);
    if(signarg[i]==0){ skewterm[i] = 1; }
    if(signarg[i] > 0){ skewterm[i] = *skew2; }
    if(signarg[i] < 0){ skewterm[i] = pow(*skew2,-1); }
    u[i] = *dfpluss1 * signarg[i] * y[i]/(*df * exp(2 * lambda[i]) * skewterm[i] + pow(signarg[i], 2) ) - 1;
    lambdadagg[i+1] = *phi1 * lambdadagg[i] + *kappa1 * u[i] + *kappastar * signnegy[i] * (u[i]+1);
    lambda[i+1] = *omega + lambdadagg[i+1];
  }

}

extern "C" {
  void tegarchrecursion2 ( int * iN, double * omega, double * phi1, double * phi2, double * kappa1, double * kappa2, double * kappastar, double * df, double * skew2, double * dfpluss1, double * mueps, double * y, double * y2, double * signnegy, double * signarg, double * skewterm, double * lambda, double * lambda1dagg, double * lambda2dagg, double * u );
}

void tegarchrecursion2 ( int * iN, double * omega, double * phi1, double * phi2, double * kappa1, double * kappa2, double * kappastar, double * df, double * skew2, double * dfpluss1, double * mueps, double * y, double * y2, double * signnegy, double * signarg, double * skewterm, double * lambda, double * lambda1dagg, double * lambda2dagg, double * u ) {

  for (int i=0; i < *iN-1; i++) {
    signarg[i] = y[i] + *mueps * exp(lambda[i]);
    if(signarg[i]==0){ skewterm[i] = 1; }
    if(signarg[i] > 0){ skewterm[i] = *skew2; }
    if(signarg[i] < 0){ skewterm[i] = pow(*skew2,-1); }
    u[i] = *dfpluss1 * signarg[i] * y[i]/(*df * exp(2 * lambda[i]) * skewterm[i] + pow(signarg[i], 2) ) - 1;
    lambda1dagg[i+1] = *phi1 * lambda1dagg[i] + *kappa1 * u[i];
    lambda2dagg[i+1] = *phi2 * lambda2dagg[i] + *kappa2 * u[i] + *kappastar * signnegy[i]*(u[i]+1);
    lambda[i+1] = *omega + lambda1dagg[i+1] + lambda2dagg[i+1];
  }

}
