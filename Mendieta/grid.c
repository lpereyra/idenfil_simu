#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "grid.h"

extern void grid_init(void)
{
  unsigned long nalloc = (grid.ngrid*grid.ngrid*grid.ngrid) + 1;

  printf("allocating %.5f gb\n",
                     (double)(nalloc*sizeof(int))/1024.0/1024.0/1024.0);
  grid.icell	= (type_int *) malloc(nalloc*sizeof(type_int));
  assert(grid.icell != NULL);
  for(unsigned long i = 0; i<=nalloc; i++)
    grid.icell[i] = cp.npart;

}

extern void grid_build(void)
{
  const long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;
  long i, j, ix, iy, iz, ibox;
  type_int *Id;
  double fac;

  fac = (double)grid.ngrid/(double)cp.lbox ;
  fprintf(stdout,"Building Grid..... Ngrid = %lu\n",grid.ngrid);
  Id = (type_int *) malloc(grid.nobj*sizeof(type_int));
  assert(Id != NULL);

  for(i = 0; i < grid.nobj; i++)
  {
#ifdef COLUMN
    ix = (long)((double)P.x[i]*fac);
    iy = (long)((double)P.y[i]*fac);
    iz = (long)((double)P.z[i]*fac);
#else
    ix = (long)((double)P[i].pos[0]*fac);
    iy = (long)((double)P[i].pos[1]*fac);
    iz = (long)((double)P[i].pos[2]*fac);
#endif

    #ifdef PERIODIC
      ibox = igrid(( (ix >= (long)grid.ngrid) ? ix-(long)grid.ngrid : ( (ix<0) ? ix + (long)grid.ngrid : ix ) ),\
                   ( (iy >= (long)grid.ngrid) ? iy-(long)grid.ngrid : ( (iy<0) ? iy + (long)grid.ngrid : iy ) ),\
                   ( (iz >= (long)grid.ngrid) ? iz-(long)grid.ngrid : ( (iz<0) ? iz + (long)grid.ngrid : iz ) ),\
                   (long)grid.ngrid);
    #else
      ibox = igrid(( (ix >= (long)grid.ngrid) ? (long)grid.ngrid-1 : ( (ix<0) ? 0 : ix ) ),\
                   ( (iy >= (long)grid.ngrid) ? (long)grid.ngrid-1 : ( (iy<0) ? 0 : iy ) ),\
                   ( (iz >= (long)grid.ngrid) ? (long)grid.ngrid-1 : ( (iz<0) ? 0 : iz ) ),\
                   (long)grid.ngrid);
    #endif

    assert(ibox<nalloc);

    Id[i] = grid.icell[ibox];
    grid.icell[ibox] = i;
  }

  j = 0;
  for(ibox = 0; ibox < nalloc; ibox++)
  {      
    i = grid.icell[ibox];
    grid.icell[ibox] = j;

    while(i != cp.npart)
    {
      unsigned long tmp = Id[i];
      Id[i] = j;
      j++;	
      i = tmp;
    }
  }

  //Ghost cell
  grid.icell[nalloc] = cp.npart;

  for(i=0;i<cp.npart;i++)
  {
    if(Id[i] != i)
    {

#ifdef COLUMN

      type_real Px_source  = P.x[i];
      type_real Py_source  = P.y[i];
      type_real Pz_source  = P.z[i];
      #ifdef STORE_VELOCITIES
        type_real Pvx_source = P.vx[i];
        type_real Pvy_source = P.vy[i];
        type_real Pvz_source = P.vz[i];
      #endif
      #ifdef STORE_IDS
        type_int Pid_source  = P.id[i];
      #endif
      type_int Psub_source  = P.sub[i];

#else

      struct particle_data P_source = P[i];

#endif

      type_int  idsource = Id[i];
      type_int  dest     = Id[i];

      while(1)
      {
#ifdef COLUMN
        type_real Px_save = P.x[dest];
        type_real Py_save = P.y[dest];
        type_real Pz_save = P.z[dest];
        #ifdef STORE_VELOCITIES
          type_real Pvx_save = P.vx[dest];
          type_real Pvy_save = P.vy[dest];
          type_real Pvz_save = P.vz[dest];
        #endif
        #ifdef STORE_IDS
          type_int  Pid_save = P.id[dest];
        #endif
        type_int  Psub_save = P.sub[dest];
#else
        struct particle_data P_save = P[dest];
#endif

	      type_int idsave = Id[dest];

#ifdef COLUMN
        P.x[dest]  = Px_source;
        P.y[dest]  = Py_source;
        P.z[dest]  = Pz_source;
        #ifdef STORE_VELOCITIES
          P.vx[dest] = Pvx_source;
          P.vy[dest] = Pvy_source;
          P.vz[dest] = Pvz_source;
        #endif
        #ifdef STORE_IDS
          P.id[dest] = Pid_source;
        #endif
        P.sub[dest] = Psub_source;
#else
        P[dest] = P_source;
#endif

        Id[dest] = idsource;

        if(dest == i)  break;

#ifdef COLUMN
        Px_source  = Px_save;
        Py_source  = Py_save;
        Pz_source  = Pz_save;
  	#ifdef STORE_VELOCITIES
          Pvx_source = Pvx_save;
          Pvy_source = Pvy_save;
          Pvz_source = Pvz_save;
	#endif
	#ifdef STORE_IDS
          Pid_source = Pid_save;
	#endif
        Psub_source = Psub_save;

#else
        P_source = P_save;
#endif

        idsource = idsave;
  	    dest = idsource;
      } // cierra el while
    }  // cierra el if

  } // cierra el for


  free(Id);
}

extern void grid_free(void)
{
  if(grid.icell!=NULL) free(grid.icell);
}
