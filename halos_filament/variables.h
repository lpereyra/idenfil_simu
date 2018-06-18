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

struct data_list
{
  type_real  r;
  int  id;
};

struct segmentstd
{
  int   size;
  int   *list;
  int   flag;
  float razon;
  float len;
  float elong;
  float rms;

  // para la lista de halos //
  type_int nlista_halos;
  struct data_list *lista_halos;
  ///////////////////////////

} *Seg;

struct grup_data
{
  int        save;
  int        id;
  type_real  Pos[3];
  int     NumPart;

  // para la lista de filamentos //
  type_int nlista_fil;
  struct data_list *lista_fil;
  ///////////////////////////

} *Gr;

int  nfrac;
int  nbins;
int  ncil;
type_real *fof;
#ifdef FIXED_SEPARATION
type_real RLEN;
type_real RSEP;
#endif
#ifdef MCRITIC
type_real m_critica;
#endif
extern void init_variables(int argc, char **argv);

#endif
