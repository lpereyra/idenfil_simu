#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifdef CUT_IN_LEN
  #ifndef CUT_LEN
    #define CUT_LEN 10000            // en Kpc
  #endif
#endif

#ifdef SUAVIZADO
  #ifndef NITERA
    #define NITERA 10    
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

size_t size_real;
size_t size_int;

struct segmentstd
{
  type_int size;
  type_int *list;
  type_int flag;
  float Mass[2];
  float Vnodos[6];
  float razon;
  float len;
  float elong;
  float rms;
};

struct grup_data
{
  type_int   save;
  type_int   id;
  type_real  Pos[3];
  type_int   NumPart;
};

struct segmentstd *Seg;
struct segmentstd *Seg_cut;
struct grup_data  *Gr; 
struct grup_data  *Gr_cut;

extern void init_variables(int argc, char **argv);

#endif
