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


    /* set the viscosity term, add the wind term in the force calculation */
    for(j=0;j<Ndim;j++) {
      wind_visc_force(Nbody,f[j],vis,velo[j],wind[j]);
    }

    /* calculate distance from central mass */
    for(k=0;k<Nbody;k++){
      r[k] = 0.0;
    }
    for(i=0;i<Ndim;i++){
	    add_norm(Nbody,r,pos[i]);
    }
    for(k=0;k<Nbody;k++){
      r[k] = sqrt(r[k]);
    }
    // for(i=0;i<Nbody;i++){
	  //   r[i] = add_norm(Ndim,pos[i]);
    // }



    /* calculate central force */
    for(i=0;i<Nbody;i++){
      for(l=0;l<Ndim;l++){
        f[l][i] -= force(G*mass[i]*M_central,pos[l][i],r[i]);
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
    for(k=0;k<Npair;k++){
      delta_r[k] = 0.0;
    }
    for(i=0;i<Ndim;i++){
      add_norm(Npair,delta_r,delta_pos[i]);
    }
    for(k=0;k<Npair;k++){
      delta_r[k] = sqrt(delta_r[k]);
    }
    // for(i=0;i<Npair;i++){
    //   add_norm(Ndim,delta_pos[i]);
    // }






    /*
    * add pairwise forces.
    */
    k = 0;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
        Size = radius[i] + radius[j];
        collided=0;
        for(l=0;l<Ndim;l++){
          /*  flip force if close in */
          if( delta_r[k] >= Size ){
            f[l][i] -= force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            f[l][j] += force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
          } else{
            f[l][i] += force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            f[l][j] -= force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            collided=1;
          }
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