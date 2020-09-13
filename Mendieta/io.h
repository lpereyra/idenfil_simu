#ifndef IO_H
#define IO_H

#ifndef VARIABLES_H
  #include "variables.h"
#endif

#ifndef PROP_H
  #include "propiedades.h"
#endif

extern void init_variables(int argc, char **argv);
extern void write_properties(FILE *pfout, struct propiedades_st Prop);
#ifdef FILE_ASCII
  extern void write_properties_ascii(FILE *pfout, struct propiedades_st Prop);
#endif
 
#endif
