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

// Posiciones y velocidades de las partículas
struct particle_data 
{
  type_real      Pos[3];
  #ifdef STORE_VELOCITIES
  type_real      Vel[3];
  #endif
  type_int       sub;
};

struct grup_data
{
  type_int   save;
  type_int   id;
  type_real  Pos[3];
  type_int   NumPart;
};

extern size_t size_real;
extern size_t size_int;
extern type_int  nfrac;
extern type_real  *fof;
extern struct particle_data *P;
extern struct grup_data *Gr;

#ifdef MCRITIC
  extern type_real m_critica;
#endif

extern void init_variables(int argc, char **argv);

#endif
