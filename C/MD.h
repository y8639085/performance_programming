#ifndef MD_H
#define MD_H

#include <math.h>

inline void wind_visc_force(int N, double * restrict f, double * restrict vis, double * restrict velo, double wind) {
  int i;
  #pragma simd
  for (i=0;i<N;i++) {
    f[i] = -vis[i] * velo[i] -vis[i] * wind;
  }
}

inline double force(double W, double delta, double r){
  return W*delta/(pow(r,3.0));
}

#endif