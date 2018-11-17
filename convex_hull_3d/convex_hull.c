#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "timer.h"
#include "colores.h"
#include "leesnap.h"

int main(int argc, char **argv)
{
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  //Lee archivos de la simulacion
  read_gadget();

  select_particles_fof(fof[0]);

  read_grup_fof(fof[1]);

  free(P);
  free(Gr);

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}
