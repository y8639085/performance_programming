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

DEF double wind[Ndim];
DEF double pos[Ndim][Nbody], velo[Ndim][Nbody];
DEF double f[Ndim][Nbody];
DEF double delta_pos[3][Nbody*Nbody];
DEF double delta_r[Nbody*Nbody];
DEF double r[Nbody], vis[Nbody], mass[Nbody], radius[Nbody];
DEF int collisions;

#define G 2.0
#define M_central 1000.0
#define GxM_central 2000.0
#define PADDING 64

void evolve(int Nstep, double dt);
