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

int main(int argc, char **argv)
{
  type_int  NNN;
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  leeheader();

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  fprintf(stdout,"lbox %g Kpc\n",cp.lbox);
  GREEN("**********************************\n");

  #ifdef ORIGINAL
  read_grup_fof(fof);
  #endif
  read_segment(NNN,fof);

  free(Seg);
  free(P);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}
