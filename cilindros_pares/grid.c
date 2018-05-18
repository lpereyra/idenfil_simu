#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.h"
#include "grid.h"
#include "cosmoparam.h"

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

  for( i = 0 ; i < grid.nobj ; i++ )
  {

    ix = (unsigned long)((double)P[i].Pos[0]*fac);
    iy = (unsigned long)((double)P[i].Pos[1]*fac);
    iz = (unsigned long)((double)P[i].Pos[2]*fac);

		ibox = (ix * grid.ngrid + iy) * grid.ngrid + iz;

    grid.ll[i] = grid.llirst[ibox];
    grid.llirst[ibox] = i;
  }

#ifdef REORDER

  int j = 0;
  for(ibox=0;ibox<grid.ngrid*grid.ngrid*grid.ngrid;++ibox)
  {
    i = grid.llirst[ibox];
    grid.llirst[ibox] = -1;

    while(i != -1)
    {
      int tmp = grid.ll[i];
      grid.ll[i] = j;
      i = tmp;
      j++;
    }
  }

  struct particle_data Psave, Psource;
  unsigned int idsource, idsave, dest;

  for(i=0;i<cp.npart;i++)
  {
    if(grid.ll[i] != i)
	  {
      Psource  = P[i];
	    idsource = grid.ll[i];
      dest     = grid.ll[i];

      while(1)
      {
	       Psave  = P[dest];
	       idsave = grid.ll[dest];

	       P[dest] = Psource;
	       grid.ll[dest] = idsource;

	       if(dest == i)  break;

	       Psource = Psave;
	       idsource = idsave;
	       dest = idsource;
	    }
	  }

    ix = (unsigned long)((double)P[i].Pos[0]*fac);
    iy = (unsigned long)((double)P[i].Pos[1]*fac);
    iz = (unsigned long)((double)P[i].Pos[2]*fac);
		ibox = (ix * grid.ngrid + iy) * grid.ngrid + iz;

    grid.ll[i] = grid.llirst[ibox];
    grid.llirst[ibox] = i;
  }

#endif

}

void grid_free(void)
{
  if(grid.ll!=NULL)free(grid.ll);
  if(grid.llirst!=NULL)free(grid.llirst);
}
