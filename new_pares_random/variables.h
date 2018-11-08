#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifdef CUT_IN_LEN

  #ifndef LEN_MIN
    #define LEN_MIN 9000
  #endif
  
  #ifndef LEN_MAX
    #define LEN_MAX 11000
  #endif

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

extern size_t size_real;
extern size_t size_int;

struct segmentstd
{
  type_int size;
  type_int flag;
  type_int list[2];
  type_real Mass[2];
  type_real Vnodos[6];
  type_real razon;
  type_real len;
  type_real elong;
  type_real rms;
};

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
};

extern struct segmentstd *Seg;
extern struct grup_data *Gr;

extern type_int  nfrac;
extern type_real *fof;
#ifdef MCRITIC
extern type_real m_critica;
#endif

extern void init_variables(int argc, char **argv);

#endif
