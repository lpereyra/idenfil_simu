#ifndef CALCULA_H
#define CALCULA_H

void propiedades_halos(type_real *fof);
void calc_part(type_real xc, type_real yc, type_real zc, int icent, int *npart, type_real **mean, type_real **quad, \
               type_real *mean_rad, type_real *quad_rad, type_real *rcil2, int binsper);

#ifdef CALCULA_MEDIA
  void calc_media(type_real xc, type_real yc, type_real zc, int *numpart, type_real *vel_media, type_real *rsph2, int nsph);
#endif

int point_inside(type_real dot, type_real rlong_2);

#ifdef LOG_BINS
  void logspace(type_real *rcil2, type_real max, type_real min, int bins); 
#else
  void linspace(type_real *rcil2, type_real max, type_real min, int bins); 
#endif

#endif
