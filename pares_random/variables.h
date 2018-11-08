#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef LEN_MIN
  #define LEN_MIN 9000
#endif

#ifndef LEN_MAX
  #define LEN_MAX 11000
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

struct grup_data
{
  type_int   save;
  type_int   id;
  type_real  Pos[3];
  type_int   NumPart;
  type_real  Mass;
  type_real  aa;
  type_real  bb;
  type_real  cc;
  type_real  ParP;
  type_real  vcm[3];
  type_real  L[3];
  type_real  evec[9];
} *Gr;

int  nfrac;
type_real *fof;
type_real pmin[3], pmax[3];
#ifdef MCRITIC
type_real m_critica;
#endif
 
void init_variables(int argc, char **argv);

#endif
