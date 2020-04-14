/*
 *  Simple molecular dynamics code.
 *
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"
#include "MD.h"

void evolve(int count,double dt){
int  step;
int i,j,k,l;
int collided;
double Size;
double r, delta_r;
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








    /* calculate pairwise separation of particles */
    k = 0;
    for(i=0;i<Nbody;i++){
      #pragma simd
      for(j=i+1;j<Nbody;j++){
        for(l=0;l<Ndim;l++){
          delta_pos[k][l] = pos[i][l] - pos[j][l];
        }
        k = k + 1;
      }
    }

    /*
    * add pairwise forces.
    */
    k = 0;
    double G_ij;
    double force_result;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
        /* calculate norm of separation vector */
        delta_r = add_norm(Ndim, delta_pos[k]);
        Size = radius[i] + radius[j];
        collided=0;
        /*  flip force if close in */
        G_ij = G*mass[i]*mass[j];
        if( delta_r >= Size ){
          #pragma simd
          for(l=0;l<Ndim;l++){
            f[i][l] -= force(G_ij,delta_pos[k][l],delta_r);
            f[j][l] += force(G_ij,delta_pos[k][l],delta_r);
          }
        }
        else{
          #pragma simd
          for(l=0;l<Ndim;l++){
            f[i][l] += force(G_ij,delta_pos[k][l],delta_r);
            f[j][l] -= force(G_ij,delta_pos[k][l],delta_r);
          }
          collided=1;
        }
        if( collided == 1 ){
          collisions++;
        }
        k = k + 1;
      }
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