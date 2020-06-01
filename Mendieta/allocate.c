#include <stdlib.h>
#include <stdio.h>
#include "variables.h"
#include "allocate.h"

#ifdef COLUMN
  extern int allocate_particles(struct particle_data  *Q, const type_int size)
#else
  extern int allocate_particles(struct particle_data **Q, const type_int size)
#endif
{

#ifdef COLUMN
  Q->x  = NULL;
  Q->y  = NULL;
  Q->z  = NULL;

  Q->x  = (type_real *) malloc(size*sizeof(type_real));
  Q->y  = (type_real *) malloc(size*sizeof(type_real));
  Q->z  = (type_real *) malloc(size*sizeof(type_real));

  if(!Q->x || !Q->y || !Q->z) 
  {
    fprintf(stderr, "cannot allocate pos particles\n" );
    return(0);
  }

  #ifdef STORE_VELOCITIES
  Q->vx = NULL;
  Q->vy = NULL;
  Q->vz = NULL;

  Q->vx = (type_real *) malloc(size*sizeof(type_real));
  Q->vy = (type_real *) malloc(size*sizeof(type_real));
  Q->vz = (type_real *) malloc(size*sizeof(type_real));

  if(!Q->vx || !Q->vy || !Q->vz) 
  {
    fprintf(stderr, "cannot allocate vel particles\n" );
    return(0);
  }
  #endif

  #ifdef STORE_IDS
  Q->id = NULL;
  Q->id = (type_int  *) malloc(size*sizeof(type_int));
  
  if(!Q->id) 
  {
    fprintf(stderr, "cannot allocate ids particles\n" );
    return(0);
  }
  #endif

#else

  *Q = NULL;
  *Q = (struct particle_data *) malloc(size*sizeof(struct particle_data));

  if(!*Q) 
  {
    fprintf(stderr, "cannot allocate particles\n" );
    return(0);
  }    

#endif

  return ( 1 );
}

#ifdef COLUMN
  extern void free_particles(struct particle_data  *Q)
#else
  extern void free_particles(struct particle_data **Q)
#endif
{

#ifdef COLUMN
  if(Q->x)  free(Q->x);
  if(Q->y)  free(Q->y);
  if(Q->z)  free(Q->z);
  #ifdef STORE_VELOCITIES
  if(Q->vx) free(Q->vx);
  if(Q->vy) free(Q->vy);
  if(Q->vz) free(Q->vz);
  #endif   
  #ifdef STORE_IDS
  if(Q->id) free(Q->id);
  #endif   
#else
  if(*Q) free(*Q);
#endif  

}

