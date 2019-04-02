#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifdef ITERA
  #ifndef NITERA
    #define NITERA 1
  #endif
#endif

#ifdef FIX_NSMOOTH
  #ifndef NSMOOTH
    #define NSMOOTH 30
  #endif
#else
  #ifndef RSPACE
    #define RSPACE 500.0
  #endif
#endif

#ifndef RVIR_FACTOR
  #define RVIR_FACTOR 2.0
#endif

#ifndef R_SPH
  #define R_SPH 2000.0
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

#ifdef PARTICLES
  struct particle_data 
  {
    type_real      Pos[3];
    type_real      Vel[3];
    #ifdef STORE_IDS
    type_int       id;
    #endif
  } *P;
#endif

struct segmentstd
{
  type_int   size;
  type_real *Pos_list;
  type_int flag;
  type_real Mass[2];
  type_real Vnodos[6];
  //#ifdef CALCULA_VCM
    type_real Vmedia[3];
  //#endif
  type_real razon;
  type_real len;
  type_real elong;
  type_real rms;
  #ifdef PARTICLES
  type_real vol;
  type_real rho;
  type_real mu;
  type_real sigma_per;
  type_real mass_part;
  type_real Rvir[2];
  #endif
} *Seg;

struct grup_data
{
  type_int   save;
  type_int   id;
  type_real  Pos[3];
  type_int   NumPart;
} *Gr;

type_int  nfrac;
type_real *fof;
#ifdef MCRITIC
type_real m_critica;
#endif

extern void init_variables(int argc, char **argv);

#endif
