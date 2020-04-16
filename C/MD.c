/*
 *  Simple molecular dynamics code.
 *
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"
#include "MD.h"

void evolve(int count,double dt){
int step;
int i,j,k,l;
double Size;
double r, delta_r;
double force_val;
double f2[Ndim + PADDING];
  /*
  * Loop over timesteps.
  */
  for(step = 1;step<=count;step++){
    printf("timestep %d\n",step);
    printf("collisions %d\n",collisions);

    /* set the viscosity and the wind term in the force calculation */
    #pragma simd
    for(j=0;j<Nbody;j++) {
      wind_visc_force(Ndim,f[j],vis[j],velo[j],wind);

    /* calculate distance from central mass */
      r = add_norm(Ndim, pos[j]);
      /* calculate central force */
      for(l=0;l<Ndim;l++){
        f[j][l] -= force(GxM_central*mass[j],pos[j][l],r);
      }
    }

    /*
    * add pairwise forces.
    */
    k = 0;
    double G_ij;
    for(i=0;i<Nbody;i++){
      memset(f2, 0, Ndim*sizeof(double));
      for(j=i+1;j<Nbody;j++){
        /* calculate pairwise separation of particles */
        for(l=0;l<Ndim;l++){
          delta_pos[l] = pos[i][l] - pos[j][l];
        }
        /* calculate norm of separation vector */
        delta_r = add_norm(Ndim, delta_pos);
        Size = radius[i] + radius[j];
        /*  flip force if close in */
        G_ij = G*mass[i]*mass[j];
        if( delta_r >= Size ){
          #pragma simd
          for(l=0;l<Ndim;l++){
            // force_val = force(G_ij,delta_pos[l],delta_r);
            f2[l] -= force(G_ij,delta_pos[l],delta_r);
            f[j][l] += force(G_ij,delta_pos[l],delta_r);
          }
        }
        else{
          #pragma simd
          for(l=0;l<Ndim;l++){
            // force_val = force(G_ij,delta_pos[l],delta_r);
            f2[l] += force(G_ij,delta_pos[l],delta_r);;
            f[j][l] -= force(G_ij,delta_pos[l],delta_r);;
          }
          collisions++;
        }
        k = k + 1;
      }
      for(l=0; l<Ndim; l++)
        f[i][l] += f2[l];
    }

    /* update positions and velocities */
    #pragma simd
    for(i=0;i<Nbody;i++){
      for(j=0;j<Ndim;j++){
        pos[i][j] += dt * velo[i][j];
        velo[i][j] += dt * (f[i][j]/mass[i]);
      }
    }
  }
}