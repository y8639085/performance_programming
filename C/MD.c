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
  /*
  * Loop over timesteps.
  */
  for(step = 1;step<=count;step++){
    printf("timestep %d\n",step);
    printf("collisions %d\n",collisions);

    /* set the viscosity and the wind term in the force calculation */
    for(j=0;j<Ndim;j++) {
      wind_visc_force(Nbody,f[j],vis,velo[j],wind[j]);
    }

    /* calculate distance from central mass */
    memset(r, 0.0, Nbody*sizeof(double));
    for(k=0;k<Nbody;k++){ 
      for(i=0;i<Ndim;i++){
        r[k] += (pos[i][k]* pos[i][k]);// inline and vectorisation
      }
      r[k] = sqrt(r[k]);
    }




    /* calculate central force */
    for(i=0;i<Nbody;i++){
      for(l=0;l<Ndim;l++){
        f[l][i] -= force(GxM_central*mass[i],pos[l][i],r[i]);
      }
    }
    /* calculate pairwise separation of particles */
    k = 0;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
        for(l=0;l<Ndim;l++){
          delta_pos[l][k] = pos[l][i] - pos[l][j];
        }
        k = k + 1;
      }
    }

    /* calculate norm of separation vector */
    memset (delta_r, 0.0, Npair * sizeof (double));
    for(k=0;k<Npair;k++){
      for (i = 0; i < Ndim; i++) {
        delta_r[k] += (delta_pos[i][k] * delta_pos[i][k]);
      }
      delta_r[k] = sqrt(delta_r[k]);
    }





    /*
    * add pairwise forces.
    */
    k = 0;
    double G_ij;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
        Size = radius[i] + radius[j];
        collided=0;
        /*  flip force if close in */
        G_ij = G*mass[i]*mass[j];
        if( delta_r[k] >= Size ){
          for(l=0;l<Ndim;l++){
            f[l][i] -= force(G_ij,delta_pos[l][k],delta_r[k]);
            f[l][j] += force(G_ij,delta_pos[l][k],delta_r[k]);
          }
        }
        else{
          for(l=0;l<Ndim;l++){
            f[l][i] += force(G_ij,delta_pos[l][k],delta_r[k]);
            f[l][j] -= force(G_ij,delta_pos[l][k],delta_r[k]);
          }
          collided=1;
        }
        if( collided == 1 ){
          collisions++;
        }
        k = k + 1;
      }
    }

    /* update positions */
    for(i=0;i<Nbody;i++){
      for(j=0;j<Ndim;j++){
        pos[j][i] += dt * velo[j][i];
      }
    }

    /* update velocities */
    for(i=0;i<Nbody;i++){
      for(j=0;j<Ndim;j++){
        velo[j][i] += dt * (f[j][i]/mass[i]);
      }
    }
  }
}