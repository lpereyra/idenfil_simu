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

struct grup_data
{
  int        save;
  type_int   id;
  type_real  Pos[3];
  int        NumPart;
  int        *list;
  #ifdef CALCULA_MEDIA
  type_real   *Vmedia;
  #endif
} *Gr;

int  nfrac;
int  nbins;
type_real *fof;
type_real pmin[3], pmax[3];
type_real len_min, len_max;
#ifdef MCRITIC
type_real m_critica;
#endif
 
void init_variables(int argc, char **argv);

#endif
