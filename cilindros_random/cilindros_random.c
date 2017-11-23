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
#include "calcula_random.h"

int main(int argc, char **argv)
{
  int    i,NNN;
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);
  // Lee archivos de la simulacion 
  read_gadget();

  GREEN("********** Important *************\n");
  sprintf(message,"Mpart %g\n",cp.Mpart*1.e10);RED(message);
  sprintf(message,"Nodos Mass %g Npart %d\n",NNN*cp.Mpart*1e10,NNN);RED(message);
  GREEN("**********************************\n");
  fflush(stdout);

  read_segment_random(NNN,fof);

  /////////////////////////////////////////////////////////////////////////
  propiedades_random(NNN,fof);
  /////////////////////////////////////////////////////////////////////////
  free(P);

  for(i=0;i<cp.nseg;i++)
  {
    free(Segrand[i].Pos);
    free(Segrand[i].Vmedia);
  }
  free(Segrand);
  free(PropSeg);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}
