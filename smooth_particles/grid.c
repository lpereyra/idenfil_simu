#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "grid.h"
#include "cosmoparam.h"

extern void grid_init()
{
  unsigned long i;
  const unsigned long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;

  printf("allocating %.5f gb\n",
                     (double)((grid.nobj+2*nalloc)*sizeof(int))/1024.0/1024.0/1024.0);
  grid.icell	= (type_int *) malloc(nalloc*sizeof(type_int));
  assert(grid.icell != NULL);
  grid.size	= (type_int *) malloc(nalloc*sizeof(type_int));
  assert(grid.size != NULL);

  grid.ll = (type_int *) malloc(grid.nobj*sizeof(type_int));
  assert(grid.ll != NULL);

  for(i = 0; i<nalloc; i++)
  {
    grid.icell[i] = grid.nobj;
    grid.size[i] = 0;
  }
  return;
}

inline unsigned long igrid(const long ix, const long iy, const long iz, const long ngrid)
{
  return (ix * ngrid + iy) * ngrid + iz;
}

extern void grid_build(type_int *id_fil, type_int *k_fil)
{
  const long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;
  long i, ix, iy, iz, ibox;
  type_int idx, kk;
  double fac;

  fac = (double)grid.ngrid/(double)cp.lbox ;
  fprintf(stdout,"Building Grid..... Ngrid = %lu\n",grid.ngrid);

  i = 0;
  for(idx = 0; idx < cp.nseg; idx++)
  {
    for(kk=0;kk<Seg[idx].size;kk++)
    {
      ix = (long)((double)Seg[idx].Pos_list[3*kk+0]*fac);
      iy = (long)((double)Seg[idx].Pos_list[3*kk+1]*fac);
      iz = (long)((double)Seg[idx].Pos_list[3*kk+2]*fac);

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

      grid.ll[i] = grid.icell[ibox];
      grid.icell[ibox]  = i;
      grid.size[ibox]  += 1;

      id_fil[i] = idx;
      k_fil[i]  = kk;
      i++;
    }
  }

  return;

}

extern void grid_free(void)
{
  if(grid.icell!=NULL) free(grid.icell);
  if(grid.ll!=NULL)    free(grid.ll);

  return;
}
