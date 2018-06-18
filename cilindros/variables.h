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
#ifdef COLUMN

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

#else

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

#endif

struct segmentstd
{
  type_int  size;
  type_int  * list;
  type_int  flag;
  type_real razon;
  type_real len;
  type_real elong;
  type_real rms;
  type_real *Vmedia;
} *Seg;

struct grup_data
{
  type_int  save;
  type_int  id;
  type_real Pos[3];
  #ifdef VEL_RELATIVA
  type_real  Vnod[3];
  #endif
  type_int  NumPart;
  type_int  * list;
} *Gr;

type_int  nfrac;
type_int  nbins;
type_int  ncil;
type_real *fof;
#ifdef FIXED_SEPARATION
type_real RLEN;
type_real RSEP;
#endif
#ifdef MCRITIC
type_real m_critica;
#endif

extern void init_variables(int argc, char **argv);
#ifdef COLUMN
  extern int allocate_particles(const type_int size);
  extern void free_particles( void );
#endif

#endif
