#include <stdlib.h>
#include <stdio.h>
#include "variables.h"
#include "allocate.h"

/*  file allocate.c
 *  Routines for allocating particle
*/

/*
 *  Allocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
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

  Q->sub = NULL;

  Q->sub = (type_int  *) malloc(size*sizeof(type_int));
  
  if(!Q->sub) 
  {
    fprintf(stderr, "cannot allocate sub particles\n" );
    return(0);
  }

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

/*
 *  Deallocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
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
  if(Q->sub) free(Q->sub);
#else
  if(*Q) free(*Q);
#endif  

}

/*
 *  Reallocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
#ifdef COLUMN
  extern int reallocate_particles(struct particle_data  *Q, const type_int size)
#else
  extern int reallocate_particles(struct particle_data **Q, const type_int size)
#endif
{

#ifdef COLUMN
  Q->x  = (type_real *) realloc(Q->x,size*sizeof(type_real));
  Q->y  = (type_real *) realloc(Q->y,size*sizeof(type_real));
  Q->z  = (type_real *) realloc(Q->z,size*sizeof(type_real));

  if(!Q->x || !Q->y || !Q->z) 
  {
    fprintf(stderr, "cannot reallocate pos particles\n" );
    return(0);
  }

  #ifdef STORE_VELOCITIES
  Q->vx = (type_real *) realloc(Q->vx,size*sizeof(type_real));
  Q->vy = (type_real *) realloc(Q->vy,size*sizeof(type_real));
  Q->vz = (type_real *) realloc(Q->vz,size*sizeof(type_real));

  if(!Q->vx || !Q->vy || !Q->vz) 
  {
    fprintf(stderr, "cannot reallocate vel particles\n" );
    return(0);
  }
  #endif

  #ifdef STORE_IDS
  Q->id = (type_int  *) realloc(Q->id,size*sizeof(type_int));
  
  if(!Q->id) 
  {
    fprintf(stderr, "cannot reallocate ids particles\n" );
    return(0);
  }
  #endif

  Q->sub = (type_int  *) realloc(Q->sub,size*sizeof(type_int));
  
  if(!Q->sub) 
  {
    fprintf(stderr, "cannot reallocate sub particles\n" );
    return(0);
  }

#else

  *Q = (struct particle_data *) realloc(*Q,size*sizeof(struct particle_data));

  if(!*Q) 
  {
    fprintf(stderr, "cannot reallocate particles\n" );
    return(0);
  }    

#endif

  return ( 1 );
}


