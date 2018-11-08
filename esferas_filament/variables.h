#ifndef VARIABLES_H
#define VARIABLES_H

#include <vector>
#include <set>
#include <utility>

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef RAUX
  #define RAUX 5000.0                  // en Kpc
#endif
 
#ifdef SPH

  #ifndef H_SOFT
    #define H_SOFT 1000.0                  // en Kpc
  #endif

  #define H_SOFT_2 H_SOFT*H_SOFT           // en Kpc

#else
 
  #ifndef K_NEIGH_SIZE
    #define K_NEIGH_SIZE 32
  #endif

  #define SET_MAX_SIZE 2*K_NEIGH_SIZE
 
#endif  

#define N_part_types 6    /* Number of particle types */

/* Precision del codigo (reales) */
#ifdef PRECDOUBLE
typedef double type_real;
#else
typedef float type_real;
#endif

/* Precision del codigo (enteros) */
#ifdef LONGIDS
typedef unsigned long long type_int;
#else
typedef unsigned int type_int;
#endif

extern size_t size_real;
extern size_t size_int;

#ifdef PARTICLES
  struct particle_data 
  {
    type_real      Pos[3];
    #ifdef STORE_VELOCITIES
    type_real      Vel[3];
    #endif
    #ifdef STORE_IDS
    type_int       id;
    #endif
  };
#endif


struct segmentstd
{
  type_int size;
  type_int start;
  type_int flag;
  type_real Mass[2];
  type_real Vnodos[6];
  type_real razon;
  type_real len;
  type_real elong;
  type_real rms;
  type_real *dens;
};

struct grup_data
{
  type_int   save;
  type_int   id;
  type_real  Pos[3];
  type_int   NumPart;
};

struct Seg_centr_data
{
  type_int   save;
  type_int   id;
  type_real  Pos[3];
};


#ifdef PARTICLES
extern struct particle_data *P;
#endif
extern struct segmentstd *Seg;
extern struct grup_data *Gr;
extern struct Seg_centr_data *Seg_centr;

extern type_int  nfrac;
extern type_real *fof;
#ifdef MCRITIC
extern type_real m_critica;
#endif

extern void init_variables(int argc, char **argv);

#endif
