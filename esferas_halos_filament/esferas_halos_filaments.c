#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "colores.h"
#include "calcula.h"
#include "grid.h"

struct gridst grid;

int main(int argc, char **argv)
{
  type_int NNN;
  double start,end;
  const type_real r_aux = 2.0f*RAUX; // RADIO DE BUSQUEDA EN Kpc

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  #ifdef PARTICLES
  read_gadget();
  #else
  leeheader();
  #endif
  
  #if defined(ORIGINAL) || !defined(PARTICLES)
    read_grup_fof(fof);
  #endif
  read_segment(NNN,fof);

  #if defined(ORIGINAL) && defined(PARTICLES)
  free(Gr);
  #endif

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  fprintf(stdout,"lbox %g Kpc\n",cp.lbox);
  //sprintf(filename,"lbox %g Kpc\n",cp.lbox);GREEN(filename);
  GREEN("**********************************\n");

  grid.nobj = cp.ncent_seg;
  grid.ngrid = (long)(cp.lbox/r_aux);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }else{
    fprintf(stdout,"Using NGRID = %lu\n",grid.ngrid);
  }

  fflush(stdout);

  grid_init();
  grid_build();

  propiedades(NNN,fof,r_aux);
  /////////////////////////////////////////////////////////////////////////

  grid_free();
  free(Seg);
  free(Seg_centr);
  #if defined(ORIGINAL) || !defined(PARTICLES)
  free(Gr);
  #endif
  #if defined(PARTICLES)
  free(P);
  #endif 

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}
