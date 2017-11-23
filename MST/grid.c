#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.h"
#include "grid.h"
#include "cosmoparam.h"

struct gridst grid;

void grid_init(void)
{
  unsigned long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;

  printf("allocating %.5f gb\n",
                     (double)((grid.nobj+nalloc)*sizeof(int))/1024.0/1024.0/1024.0);
	grid.llirst	= (int *) malloc(nalloc*sizeof(int));
  assert(grid.llirst != NULL);
  memset(grid.llirst,-1,nalloc*sizeof(int));
	grid.ll = (int *) malloc(grid.nobj*sizeof(int));
  assert(grid.ll != NULL);
}

void grid_build(void)
{
  int i;
  unsigned long ix, iy, iz;
  double fac;
	unsigned long ibox;

  fac = (double)grid.ngrid/(double)cp.lbox ;
	printf("Building Grid..... Ngrid = %lu\n",grid.ngrid);

  for( i = 0 ; i < (int)grid.nobj ; i++ )
  {

    if(grid.step!=0 && P[i].sub==0) continue;

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

void grid_free(void)
{
  if(grid.ll!=NULL)free(grid.ll);
  if(grid.llirst!=NULL)free(grid.llirst);
}
