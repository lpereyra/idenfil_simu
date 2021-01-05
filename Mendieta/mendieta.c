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

  omp_set_nested(1);

  TIMER(start);
  
  init_variables(argc,argv);

  /* Read Simulation */
  read_gadget();

#ifdef CHANGE_POSITION
  RED("\nBegins Change Positions\n");
  fflush(stdout);
  change_positions(cp.npart);
  GREEN("\nEnd Change Positions\n");
  fflush(stdout);
#endif

  fprintf(stdout, "\nBegins Identification\n");
  fflush(stdout);
  identification();
  fprintf(stdout, "\nEnds Identification\n");
  fflush(stdout);

  free_particles(&P);
  TIMER(end);
  
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}
