#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef RVIR_FACTOR
  #define RVIR_FACTOR 2.0
#endif

#ifdef CUT_IN_LEN

  #ifndef LEN_MIN
    #define LEN_MIN 10000   // en Kpc
  #endif

  #ifndef LEN_MAX
    #define LEN_MAX 100000  // en Kpc
  #endif

#endif

#ifndef TYPE_FLAG
  #define TYPE_FLAG 2
#endif

#ifdef CUT_ELONGACION
  #ifndef CUT_ELONG
    #define CUT_ELONG 0.9
  #endif
#endif

#ifndef NRANDOM
  #define NRANDOM 100
#endif

#ifndef NLIMITE
  #define NLIMITE 10*NRANDOM
#endif

#ifndef RSPH
  #define RSPH 2000
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
  type_int id;
  type_int size;
  type_int start;
  type_real Rvir_2[2];
} *Seg;

struct grup_data
{ 
  type_real  Pos[3];
} *Gr;

type_int  nfrac;
type_int  ncil;
type_real *fof;
type_real RAUX;
type_real RINIT;

extern void init_variables(int argc, char **argv);

#endif
