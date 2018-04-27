#ifndef CALCULA_H
#define CALCULA_H

#include "list.h"

void propiedades(int NNN, type_real *fof);
#ifdef SAVEPART

  void calcular_mean(const int i, const int binsper, const type_real rsep, \
           const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
  type_real **vmean, type_real **vquad, type_real *vmean_par, type_real *vquad_par,\
  type_real *vmean_perp, type_real *vquad_perp, type_real *nodo, struct node_sph **root,\
  int *size_part, int *vpart);

  void calc_part(type_real * const Pos_cent, type_real * const Vmedia, int *vec, int *ncont, int *npart, type_real **mean, type_real **quad, \
  type_real *versor, type_real *mean_par, type_real *quad_par, type_real *mean_perp, type_real *quad_perp, type_real *rsph2, type_real rlong_2, int nsph);

  #ifdef EXTEND

    void ext_mean(FILE *pfextend, const int i, const int binsper, const type_real rsep,\
    const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
    int *size_part, int *vpart);

  #endif
 
#else        

  void calcular_mean(const int i, const int binsper, const type_real rsep, \
           const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
  type_real **vmean, type_real **vquad, type_real *vmean_par, type_real *vquad_par,\
  type_real *vmean_perp, type_real *vquad_perp, type_real *nodo, struct node_sph **root);

  void calc_part(type_real * const Pos_cent, type_real * const Vmedia, int *npart, type_real **mean, type_real **quad, \
  type_real *versor, type_real *mean_par, type_real *quad_par, type_real *mean_perp, type_real *quad_perp, type_real *rsph2, type_real rlong_2, int nsph);

  #ifdef EXTEND

    void ext_mean(FILE *pfextend, const int i, const int binsper, const type_real rsep,\
    type_real rlong, type_real rlong_2, type_real * const rcil2);
 
  #endif

#endif

#ifdef CALCULA_MEDIA

  void cylinder_mean(const int i, const type_real rsep, const type_real rlong, \
                     const type_real rlong_2, type_real * const rcil2);

  void calc_media(type_real * const Pos_cent, type_real * const versor, int *numpart, type_real *vel_media, type_real * const rsph2, const type_real rlong_2, const int nsph);

#endif

int point_inside(type_real dot, type_real rlong_2);

int cmp(const void * a, const void * b);
void stadistic(int n, type_real *MAX, type_real *MIN, type_real *LMAX);
#ifdef BIN_LOG
  void logspace(type_real *rcil2, type_real max, type_real min, int bins);
#else
  void linspace(type_real *rcil2, type_real max, type_real min, int bins); 
#endif

#endif
