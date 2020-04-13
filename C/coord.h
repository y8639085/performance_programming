/*
 * This file defines static arrays that contains the primary coordinates
 * of the particles,
 *
 *  Nbody	  Number of particles
 *  Npair	  Number of particle pairs
 *  pos		  Position of the particles
 *  r         distance of partice from central mass 
 *  vel		  velocity of the particles
 *  f		  Forces acting on each particle
 *  vis       viscosity coefficient for each particle
 *  mass	  mass of each particle
 *  delta_pos separation vector for each particle pair
 *  delta_r	  separation for each particle pair
 */

#ifdef DECL
#define DEF
#else
#define DEF extern
#endif
#define Nbody 4*1024
#define  Npair ((Nbody*(Nbody-1))/2)

enum{ Xcoord=0, Ycoord, Zcoord, Ndim };

#define PADDING 64
DEF double pos[Ndim][Nbody + PADDING], velo[Ndim][Nbody + PADDING];
DEF double f[Ndim][Nbody + PADDING], vis[Nbody + PADDING], mass[Nbody + PADDING], radius[Nbody + PADDING];
DEF double delta_pos[3][Nbody*Nbody + PADDING];
DEF double r[Nbody + PADDING];
DEF double delta_r[Nbody*Nbody + PADDING];
DEF double wind[Ndim + PADDING];
DEF int collisions;

#define G 2.0
#define M_central 1000.0
#define GxM_central 2000.0


void evolve(int Nstep, double dt);