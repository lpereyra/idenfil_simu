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
#include "create_random.h"

int main(int argc, char **argv)
{
  int    i,NNN;
  double start,end;
  char filename[200];

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  GREEN("********** Important *************\n");
  sprintf(message,"Mpart %g\n",cp.Mpart*1.e10);RED(message);
  sprintf(message,"Nodos Mass %g Npart %d\n",NNN*cp.Mpart*1e10,NNN);RED(message);
  GREEN("**********************************\n");
  fflush(stdout);

  read_segment(NNN,fof);
  read_grup_fof(fof);

  /////////////////////////////////////////////////////////////////////////
  create_rand(NNN,fof);
  /////////////////////////////////////////////////////////////////////////

  free(P);

  for(i=0;i<cp.nseg;i++)
    free(Seg[i].list);
  free(Seg);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}
