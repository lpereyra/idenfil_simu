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
#include "grid.h"
#include "calcula_perfil.h"

int main(int argc, char **argv)
{
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  // Lee archivos de la simulacion 
  read_gadget();

  read_grup_fof(fof);

  /////////////////////////////////////////////////////////////////////////
  propiedades_halos(fof);
  /////////////////////////////////////////////////////////////////////////
  free(P);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}
