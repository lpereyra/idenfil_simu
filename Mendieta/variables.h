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

#ifndef FOF_OVERDENSITY
  #define FOF_OVERDENSITY 200
#endif

#ifndef NCUT
  #define NCUT 2
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

#define RHOCRIT 2.77525E11   /* Densidad crítica del Universo [Msol h² / Mpc³] */
#define GCONS 6.67300E-20    /* Constante de Gravitación [km³ / kg / seg²]     */
#define Msol 1.9891E30       /* Masa del Sol [Kg]                              */
#define Kpc 3.08568025E16    /* Kiloparsec -> Kilometro                        */  

extern struct cosmoparam
{
  double   omegam      ;  /* Omega Materia                         */
  double   omegal      ;  /* Omega Lambda                          */
  double   omegak      ;  /* Omega Curvatura                       */
  double   hparam      ;  /* Parámetro de Hubble adimensional      */
  double   lbox        ;  /* Lado del box [Kpc / h]                */
  double   Mpart       ;  /* Masa de la partícula [10^10 Msol / h] */
  type_int npart       ;  /* Número de partículas                  */
  double   redshift    ;  /* Redshift                              */
  double   aexp        ;  /*                                       */
  double   Hubble_a    ;  /*                                       */
  double   soft        ;  /* Softening [kpc / h]                   */
} cp;

/* Input and output files */
extern struct SnapST
{
  int nfiles;
  char root[200], name[200];
  int num;
} snap;

extern struct gridst
{
  unsigned long ngrid;
  unsigned long nobj;
  type_int *icell;
} grid;

/* Posiciones, velocidades y energias de las partículas */
#ifdef COLUMN

  extern struct particle_data 
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
    type_int       *sub;
  } P;

#else

  extern struct particle_data 
  {
    type_real pos[3];
    #ifdef STORE_VELOCITIES
    type_real vel[3];
    #endif
    #ifdef STORE_IDS
    type_int       id;
    #endif
    type_int       sub;
  } *P;

#endif

extern type_int  nfrac;
extern type_real *fof;
extern char message[200];
#ifdef CHANGE_POSITION
  extern type_real pmin[3], pmax[3];
#endif

#endif
