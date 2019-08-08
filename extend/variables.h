#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef RAUX
  #define RAUX 500
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

size_t size_real;
size_t size_int;

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
} *P;

struct segmentstd
{
  type_int size;
  type_int flag;
  type_real Mass[2];
  type_real Vnodos[6];
  type_real razon;
  type_real len;
  type_real elong;
  type_real rms;
  type_real *Pos;
  type_real rmin[3];
  type_real rmax[3];
  #ifdef NEW_VERSION
  type_real vol;
  type_real rho;
  type_real mu;
  type_real mass_part;
  type_real sigma;
  type_real r_pos;
  int id_pos;
  #endif
} *Seg;


#ifdef ORIGINAL

struct grup_data
{ 
  type_int   save;
  type_int   id;
  type_real  Pos[3];
  type_int   NumPart;
} *Gr;

#endif

type_int  nfrac;
type_real *fof;

extern void init_variables(int argc, char **argv);

#endif
