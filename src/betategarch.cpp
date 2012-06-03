#include <R.h>


extern "C" {
  void file1de03a27a74 ( int * iN, double * delta, double * phi1, double * kappa1, double * df, double * y2, double * dfpluss1y2, double * kappa1starsignnegy, double * lambda1, double * u );
}

void file1de03a27a74 ( int * iN, double * delta, double * phi1, double * kappa1, double * df, double * y2, double * dfpluss1y2, double * kappa1starsignnegy, double * lambda1, double * u ) {

  for (int i=0; i < *iN-1; i++) {
    u[i] = dfpluss1y2[i]/(*df * exp(2 * lambda1[i]) + y2[i]) - 1;
    lambda1[i+1] = *delta + *phi1 * lambda1[i] + *kappa1 * u[i] + kappa1starsignnegy[i] * (u[i] + 1);
  }

}

