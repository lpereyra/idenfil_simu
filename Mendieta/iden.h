#ifndef IDENTIFICADOR_H
#define IDENTIFICADOR_H

struct iden_st
{
  double *r0;        /*Linking Length para el FoF*/
  type_int nobj;
  type_int ngrupos;
};

struct temporary 
{
  type_int *head;
  type_int *ll;
  type_int *npgrup;
};

extern void identification(void);

#endif
