#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef HSPH
  #define HSPH      2000.0  // en Kpc 
#endif

#ifdef CUT_IN_LEN

  #ifndef LEN_MIN
    #define LEN_MIN 9000         // en Kpc
  #endif

  #ifndef LEN_MAX
    #define LEN_MAX 11000         // en Kpc
  #endif

#endif

#define N_part_types 6    /* Number of particle types */

// Precision del codigo (reales) //
#ifdef PRECDOUBLE
  typedef double type_real;
#else
  typedef float type_real;
#endif

// Precision del codigo (enteros) //
#ifdef LONGIDS
  typedef long type_int;
#else
  typedef unsigned int type_int;
#endif

/* Posiciones, velocidades y energias de las partículas */
struct particle_data 
{
  type_real      Pos[3];
  #ifdef STORE_VELOCITIES
  type_real      Vel[3];
  #endif
  #ifdef STORE_IDS
  type_int       id;
  #endif
  #ifndef SPH
  type_int flag;
  #endif
};

struct segmentstd
{
  type_int id;
  type_int size;
  type_int flag;
  type_real Mass[2];
  type_real Rvir[2];
  type_real Vnodos[6];
  type_real razon;
  type_real len;
  type_real elong;
  type_real rms;
  type_real *Pos;
  type_real vol;
  type_real rho;
  type_real mu;
  type_real mass_part;
};

extern type_int  nfrac;
extern type_real *fof;
extern struct particle_data *P;
extern struct segmentstd  *Seg;

extern void init_variables(int argc, char **argv);

#endif
