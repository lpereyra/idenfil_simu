#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef NX
  #define NX 2         // Numero de divisiones por coordenada.
#endif
  
#ifndef NY
  #define NY 2         // El numero de procesos tiene que ser igual
#endif

#ifndef NZ
  #define NZ 2         // al producto de ellos.
#endif

#ifndef NTHREADS
  #define NTHREADS 8   // NX*NY*NZ
#endif

#define N_part_types 6    /* Number of particle types */

// Precision del codigo (reales) //
#ifdef PRECDOUBLE
  typedef double type_real;
#else
  typedef float type_real;
#endif

// Precision del codigo (enteros) //
#ifdef LONGIDS
  typedef long type_int;
#else
  typedef int type_int;
#endif

#ifdef SAVENEIGH
  #include <vector>
#endif

/* Posiciones, velocidades y energias de las partículas */
struct particle_data 
{
  type_int       Save;
  type_real      Pos[3];
  #ifdef STORE_VELOCITIES
  type_real      Vel[3];
  #endif
  #ifdef STORE_IDS
  type_int       id;
  #endif

  double vol;
  #ifdef SAVENEIGH
  std::vector<type_int> neigh;
  #endif
  #ifdef SAVECENTROID
  double c[3];
  #endif

};

extern struct particle_data *P;

extern void init_variables(int argc, char **argv);

#endif
