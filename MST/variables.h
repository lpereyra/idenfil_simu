#ifndef VARIABLES_H
#define VARIABLES_H

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

extern int size_real;
extern int size_int;

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
  int            sub;
  int            gr;
};

struct grup_data
{
  int        save;
  type_int   id;
  type_real  Pos[3];
  int        NumPart;
};

extern int  nfrac;
extern int  ngrup;
extern type_real  *fof;
extern struct particle_data *P;
extern struct grup_data *Gr;

#ifdef BRANCH_SURVIVE
extern int N_part_survive;
#endif

#ifdef MCRITIC
extern type_real m_critica;
#endif

void init_variables(int argc, char **argv);

#endif
