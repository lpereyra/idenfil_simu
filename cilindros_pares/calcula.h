#ifndef CALCULA_H
#define CALCULA_H

void propiedades(int NNN, type_real *fof);
#ifdef SAVEPART
  void calc_part(type_real *Pos_cent, int icent, int *vec, int *ncont, int *npart, type_real **mean, type_real **quad, \
  type_real *versor, type_real *mean_par, type_real *quad_par, type_real *mean_perp, type_real *quad_perp, type_real *rsph2, type_real rlong_2, int nsph);
#else        
  void calc_part(type_real *Pos_cent, int icent, int *npart, type_real **mean, type_real **quad, \
  type_real *versor, type_real *mean_par, type_real *quad_par, type_real *mean_perp, type_real *quad_perp, type_real *rsph2, type_real rlong_2, int nsph);
#endif

#ifdef CALCULA_MEDIA
  void calc_media(type_real *Pos_cent, type_real *versor, int *numpart, type_real *vel_media, type_real *rsph2, type_real rlong_2, int nsph);
#endif

int point_inside(type_real dot, type_real rlong_2);

int cmp(const void * a, const void * b);
void stadistic(int n, type_real *MAX, type_real *MIN, type_real *LMAX); 
void logspace(type_real *rcil2, type_real max, type_real min, int bins); 
void linspace(type_real *rcil2, type_real max, type_real min, int bins); 

#endif
