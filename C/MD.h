#ifndef MD_H
#define MD_H

#include <math.h>

inline void visc_force(int N,double *f, double *vis, double *velo)
{
  int i;
  for(i=0;i<N;i++){
    f[i] = -vis[i] * velo[i];
  }
}
inline void wind_force(int N,double *f, double *vis, double velo)
{
  int i;
  for(i=0;i<N;i++){
    f[i] = f[i] -vis[i] * velo;
  }
}
inline void add_norm(int N,double *r, double *delta)
{
  int k;
  for(k=0;k<N;k++){
    r[k] += (delta[k] * delta[k]);
  }
}

// inline double add_norm(int N, double *delta)
// {
//   int k;
//   double r = 0.0;
//   for(k=0;k<N;k++){
//     r += (delta[k] * delta[k]);
//   }
//   return sqrt(r);
// }

inline double force(double W, double delta, double r){
  return W*delta/(pow(r,3.0));
}

#endif