#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "allocate.h"
#include "leesnap.h"
#include "timer.h"
#include "colores.h"
#include "iden.h"
#include "grid.h"
#include "io.h"

int main(int argc, char **argv)
{
  double start, end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  /* Read Simulation */
  read_gadget();

#ifdef CHANGE_POSITION
  RED("\nBegins Change Positions\n");
  change_positions(cp.npart);
  GREEN("\nEnd Change Positions\n");
  fflush(stdout);
#endif

  fprintf(stdout, "\nBegins Identification\n");
  identification();
  fprintf(stdout, "\nEnds Identification\n");

  /************* TERMINO LA IDENTIFICACION ***************/
  free_particles(&P);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}
