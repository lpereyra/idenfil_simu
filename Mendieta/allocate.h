#ifndef ALLOCATE_H
#define ALLOCATE_H

#ifndef VARIABLES_H
  #include "variables.h"
#endif

#ifdef COLUMN
  extern int allocate_particles(struct particle_data  *Q, const type_int size);
  extern void free_particles(struct particle_data  *Q);
#else
  extern int allocate_particles(struct particle_data **Q, const type_int size);
  extern void free_particles(struct particle_data **Q);
#endif

#endif
