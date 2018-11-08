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
    fprintf(stdout,"allocating %.5f gb\n",
                     (double)((2*nalloc)*sizeof(int))/1024.0/1024.0/1024.0);
  #else
    fprintf(stdout,"allocating %.5f gb\n",
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

  for(i=0;i<nalloc;i++)
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
    ix = (unsigned long)((double)Gr[i].Pos[0]*fac);
    iy = (unsigned long)((double)Gr[i].Pos[1]*fac);
    iz = (unsigned long)((double)Gr[i].Pos[2]*fac);

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
      
      struct grup_data Psource = Gr[i];

    	type_int  idsource = Id[i];
      type_int  dest = Id[i];

      while(1)
      {

        struct grup_data Psave = Gr[dest];

	      type_int idsave = Id[dest];

	      Gr[dest] = Psource;

        Id[dest] = idsource;

        if(dest == i)  break;

        Psource = Psave;

        idsource = idsave;
  	    dest = idsource;
      } // cierra el while
    }  // cierra el if

    ix = (unsigned long)((double)Gr[i].Pos[0]*fac);
    iy = (unsigned long)((double)Gr[i].Pos[1]*fac);
    iz = (unsigned long)((double)Gr[i].Pos[2]*fac);

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
