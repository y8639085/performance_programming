#ifndef MD_H
#define MD_H

#include <math.h>

inline void wind_visc_force(int N, double * restrict f, double vis, double * restrict velo, double * restrict wind) {
  int i;
  // #pragma ivdep
  for (i=0;i<N;i++) {
    f[i] = -vis * velo[i] -vis * wind[i];
  }
}

inline double add_norm(int N, double * restrict delta) {
  int k;
  double r = 0.0;
  for(k=0;k<N;k++){
    r += (delta[k] * delta[k]);
  }
  return sqrt(r);
}

inline double force(double W, double delta, double r){
  return W*delta/(pow(r,3.0));
}

#endif