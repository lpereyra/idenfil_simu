#ifndef CALCULA_H
#define CALCULA_H

#include "list.h"

void propiedades(int NNN, type_real *fof);
void calcular_pert(const int i, const type_real rsep, const type_real rlong, const type_real rlong_2, type_real const rcil2);
void calc_cil(const type_real * restrict Pos_cent,  const type_real * restrict versor, const type_real rcil2, const type_real rlong_2, struct list **lista);
int point_inside(type_real dot, type_real rlong_2);
void set_name(char * name, const int NNN, char * prefix);

int cmp(const void * a, const void * b);
void stadistic(int n, type_real *MAX, type_real *MIN, type_real *LMAX);
#ifdef BIN_LOG
  void logspace(type_real *rcil2, type_real max, type_real min, int bins);
#else
  void linspace(type_real *rcil2, type_real max, type_real min, int bins); 
#endif

#endif
