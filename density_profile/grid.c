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

  fprintf(stdout,"allocating %.5f gb\n",
                   (double)((grid.nobj+nalloc)*sizeof(int))/1024.0/1024.0/1024.0);

  grid.llirst	= (type_int *) malloc(nalloc*sizeof(type_int));
  assert(grid.llirst != NULL);

  grid.ll = (type_int *) malloc(grid.nobj*sizeof(type_int));
  assert(grid.ll != NULL);

  for(i = 0; i<nalloc; i++)
    grid.llirst[i] = grid.nobj;

}

inline unsigned long igrid(const long ix, const long iy, const long iz, const unsigned long ngrid)
{
  return (ix * ngrid + iy) * ngrid + iz;
}

extern void grid_build(void)
{
  unsigned long i, ix, iy, iz, ibox;
  double fac;

  fac = (double)grid.ngrid/(double)cp.lbox ;
  fprintf(stdout,"Building Grid..... Ngrid = %lu\n",grid.ngrid);

  for(i = 0; i < grid.nobj; i++)
  {
    ix = (unsigned long)((double)P[i].Pos[0]*fac);
    iy = (unsigned long)((double)P[i].Pos[1]*fac);
    iz = (unsigned long)((double)P[i].Pos[2]*fac);

    ibox = igrid(ix,iy,iz,grid.ngrid);
   
    assert(ibox>=0 && ibox<(grid.ngrid*grid.ngrid*grid.ngrid));

    grid.ll[i] = grid.llirst[ibox];
    grid.llirst[ibox] = i;
  }

}

extern void grid_free(void)
{
    if(grid.ll!=NULL)free(grid.ll);
    if(grid.llirst!=NULL)free(grid.llirst);

}
