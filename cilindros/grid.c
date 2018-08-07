#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.h"
#include "grid.h"
#include "cosmoparam.h"

extern void grid_init(void)
{
  unsigned long i;
  unsigned long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;

  #ifdef REORDER
    printf("allocating %.5f gb\n",
                     (double)((2*nalloc)*sizeof(int))/1024.0/1024.0/1024.0);
  #else
    printf("allocating %.5f gb\n",
                     (double)((grid.nobj+nalloc)*sizeof(int))/1024.0/1024.0/1024.0);
  #endif

  #ifdef REORDER
    grid.icell	= (type_int *) malloc(nalloc*sizeof(type_int));
    assert(grid.icell != NULL);

    grid.size = (type_int *) malloc(nalloc*sizeof(type_int));
    assert(grid.size != NULL);
  #else
  	grid.llirst	= (type_int *) malloc(nalloc*sizeof(type_int));
    assert(grid.llirst != NULL);

  	grid.ll = (type_int *) malloc(grid.nobj*sizeof(type_int));
    assert(grid.ll != NULL);
  #endif

  for(i = 0; i<nalloc; i++)
  {
  #ifdef REORDER
    grid.icell[i]  = grid.nobj;
    grid.size[i] = 0;
  #else
    grid.llirst[i] = grid.nobj;
  #endif
  }

}

inline unsigned long igrid(const long ix, const long iy, const long iz, const unsigned long ngrid)
{
  return (ix * ngrid + iy) * ngrid + iz;
}

extern void grid_build(void)
{
  unsigned long i, ix, iy, iz, ibox;
  double fac;
  #ifdef REORDER
  const unsigned long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;
  unsigned long j;
  type_int *Id;
  #endif

  fac = (double)grid.ngrid/(double)cp.lbox ;
  fprintf(stdout,"Building Grid..... Ngrid = %lu\n",grid.ngrid);
  #ifdef REORDER
  Id = (type_int *) malloc(grid.nobj*sizeof(type_int));
  assert(Id != NULL);
  #endif

  for(i = 0; i < grid.nobj; i++)
  {
    #ifdef COLUMN
      ix = (unsigned long)((double)P.x[i]*fac);
      iy = (unsigned long)((double)P.y[i]*fac);
      iz = (unsigned long)((double)P.z[i]*fac);
    #else
      ix = (unsigned long)((double)P[i].Pos[0]*fac);
      iy = (unsigned long)((double)P[i].Pos[1]*fac);
      iz = (unsigned long)((double)P[i].Pos[2]*fac);
    #endif

    ibox = igrid(ix,iy,iz,grid.ngrid);
      
    #ifdef REORDER
      Id[i] = grid.icell[ibox];
      grid.icell[ibox] = i;
    #else
      grid.ll[i] = grid.llirst[ibox];
      grid.llirst[ibox] = i;
    #endif
  }

#ifdef REORDER

  j = 0;
  for(ibox = 0; ibox < nalloc; ibox++)
  {      
    i = grid.icell[ibox];

    while(i != cp.npart)
    {
      unsigned long tmp = Id[i];
      Id[i] = j;
      j++;	
      i = tmp;
    }

  }

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

      #else

        struct particle_data Psource = P[i];
      

      #endif

    	type_int  idsource = Id[i];
      type_int  dest = Id[i];

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

        #else

          struct particle_data Psave = P[dest];

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

        #else

	         P[dest] = Psource;

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

	      #else

	         Psource = Psave;

	      #endif

        idsource = idsave;
  	    dest = idsource;
      } // cierra el while
    }  // cierra el if

    #ifdef COLUMN
      ix = (unsigned long)((double)P.x[i]*fac);
      iy = (unsigned long)((double)P.y[i]*fac);
      iz = (unsigned long)((double)P.z[i]*fac);
    #else
      ix = (unsigned long)((double)P[i].Pos[0]*fac);
      iy = (unsigned long)((double)P[i].Pos[1]*fac);
      iz = (unsigned long)((double)P[i].Pos[2]*fac);
    #endif

    ibox = igrid(ix,iy,iz,grid.ngrid);
    if(grid.size[ibox]==0) grid.icell[ibox] = i;
    grid.size[ibox]++;

  } // cierra el for

  free(Id);

#endif

}

extern void grid_free(void)
{
  #ifdef REORDER
    if(grid.icell!=NULL) free(grid.icell);
    if(grid.size!=NULL)  free(grid.size);
  #else
    if(grid.ll!=NULL)free(grid.ll);
    if(grid.llirst!=NULL)free(grid.llirst);
  #endif


}
