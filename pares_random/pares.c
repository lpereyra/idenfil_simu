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
#include "calcula_pares_random.h"

int main(int argc, char **argv)
{
  int NpartCut;
  double MassCut;
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  leeheader();

  read_grup_fof(fof);

  /////////////////////////////////////////////////////////////////////////

  MassCut = atof(argv[2]);
  NpartCut = (int)(MassCut/cp.Mpart) ;  // Masa de la partícula [10^10 Msol / h] 

  GREEN("********** Important *************\n");
  sprintf(message,"Mpart %g\n",cp.Mpart*1.e10);RED(message);
  sprintf(message,"Nodos Mass %g Npart %d\n",MassCut,NpartCut);RED(message);
  GREEN("**********************************\n");
  fflush(stdout);

  crea_random(NpartCut);

  /////////////////////////////////////////////////////////////////////////
  
  free(Gr);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}
