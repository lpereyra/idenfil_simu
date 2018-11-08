#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.hh"
#include "cosmoparam.hh"
#include "grid.hh"

struct gridst grid;

extern void grid_init(void)
{
  unsigned long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;

  printf("allocating %.5f gb\n",
                     (double)((grid.nobj+nalloc)*sizeof(type_int))/1024.0/1024.0/1024.0);
	grid.llirst	= (type_int *) malloc(nalloc*sizeof(type_int));
  assert(grid.llirst != NULL);
  for(unsigned long i=0;i<nalloc;i++)
    grid.llirst[i] = grid.nobj;
	grid.ll = (type_int *) malloc(grid.nobj*sizeof(type_int));
  assert(grid.ll != NULL);
}

extern void grid_build(void)
{
  unsigned long i, ix, iy, iz, ibox;
  double fac;

  fac = (double)grid.ngrid/(double)cp.lbox ;
	printf("Building Grid..... Ngrid = %lu\n",grid.ngrid);

  for( i = 0 ; i < grid.nobj ; i++ )
  {

    ix = (unsigned long)((double)P[i].Pos[0]*fac);
    iy = (unsigned long)((double)P[i].Pos[1]*fac);
    iz = (unsigned long)((double)P[i].Pos[2]*fac);

		ibox = (ix * grid.ngrid + iy) * grid.ngrid + iz;

#ifdef DEBUG
    assert(ibox >= 0L);
    if(ibox >= grid.ngrid*grid.ngrid*grid.ngrid)
      printf("%d %lu %lu %lu %lu %f %f %f\n",i,ix,iy,iz,ibox,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
    assert(ibox < grid.ngrid*grid.ngrid*grid.ngrid);
#endif

    grid.ll[i] = grid.llirst[ibox];
    grid.llirst[ibox] = i;
  }
}

extern void grid_free(void)
{
  if(grid.ll!=NULL)free(grid.ll);
  if(grid.llirst!=NULL)free(grid.llirst);
}
