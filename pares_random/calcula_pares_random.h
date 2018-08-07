#ifndef CALCULA_H
#define CALCULA_H

extern void crea_random(int NpartCut);

struct sort_prop{
  int index;
  type_real mass;
};

struct info 
{
  type_int  id_a;
  type_int  id_b;
  type_real dist;
};

struct data 
{
  type_int        nsize;
  struct info * matrix;
};

#endif
