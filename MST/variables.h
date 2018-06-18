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

extern type_int size_real;
extern type_int size_int;

/* Posiciones, velocidades y energias de las partículas */
struct particle_data 
{
  type_real      Pos[3];
  #ifdef STORE_VELOCITIES
  type_real      Vel[3];
  #endif
  //#ifdef STORE_IDS
  //type_int       id;
  //#endif
  type_int       sub;
};

struct grup_data
{
  type_int   save;
  type_int   id;
  type_real  Pos[3];
  type_int   NumPart;
};

extern type_int  nfrac;
extern type_real  *fof;
extern struct particle_data *P;
extern struct grup_data *Gr;
extern type_int *Index;

#ifdef BRANCH_SURVIVE
extern type_int N_part_survive;
#endif

#ifdef MCRITIC
extern type_real m_critica;
#endif

void init_variables(int argc, char **argv);

#endif
