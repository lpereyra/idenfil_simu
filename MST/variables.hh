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

#define RHOCRIT 2.77525E11   /* Densidad crítica del Universo [Msol h² / Mpc³] */
#define GCONS 6.67300E-20    /* Constante de Gravitación [km³ / kg / seg²]     */
#define Msol 1.9891E30       /* Masa del Sol [Kg]                              */
#define Kpc 3.08568025E16    /* Kiloparsec -> Kilometro                        */  

struct cosmoparam
{
    double  omegam			;  /* Omega Materia                         */
    double  omegal			;  /* Omega Lambda                          */
    double  omegak			;  /* Omega Curvatura                       */
		double  hparam      ;  /* Parámetro de Hubble adimensional      */
		double  lbox        ;  /* Lado del box [Kpc / h]                */
		double  Mpart       ;  /* Masa de la partícula [10^10 Msol / h] */
		type_int  npart     ;  /* Número de partículas                  */
		double  redshift		;  /* Redshift                              */
		double  aexp     		;  /*                                       */
		double  Hubble_a    ;  /*                                       */
		double  soft        ;  /* Softening [kpc / h]                   */
    type_int  ngrup     ;  /* Número de grupos                      */

};

// Posiciones y velocidades de las partículas
struct particle_data 
{
  type_real      Pos[3];
#ifdef STORE_VELOCITIES
  type_real      Vel[3];
#endif
#ifdef STORE_VELOCITIES
  type_int      Id;
#endif
  type_int       sub;
};

struct grup_data
{
  type_int   save;
  type_int   id;
  type_int   NumPart;
  type_real  Pos[3];
  type_real  aa, bb, cc;
  type_real  evec[3][3];
#ifdef STORE_VELOCITIES  
  type_real  vcm[3], sig[3], mostbound[3];
  type_real  r200, m200, v200;
  type_real  rvir, mvir, vvir;
  type_real  vmax;
  type_real  L[3]; 
  type_real  Ep, Ec, lambda;
  type_real  aa_vel, bb_vel, cc_vel;
  type_real  evec_vel[3][3];
#endif
};

extern type_int  nfrac;
extern type_real  *fof;
extern struct particle_data *P;
extern struct grup_data *Gr;
extern struct cosmoparam cp;

extern void init_variables(int argc, char **argv);

#endif
