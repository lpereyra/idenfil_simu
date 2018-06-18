#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef NPARTMIN
  #define NPARTMIN 20
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
  type_real      *x;
  type_real      *y;
  type_real      *z;
  #ifdef STORE_VELOCITIES
  type_real      *vx;
  type_real      *vy;
  type_real      *vz;
  #endif
  #ifdef STORE_IDS
  type_int       *id;
  #endif
} P;

type_int       **gr;

type_int  nfrac;
type_real *fof;
type_real pmin[3], pmax[3];

extern void init_variables(int argc, char **argv);
extern int allocate_particles(const type_int size);
extern void free_particles( void );

#endif
