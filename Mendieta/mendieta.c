#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "iden.h"
#include "colores.h"
//#include "peano.h"
#include "grid.h"

int main(int argc, char **argv)
{
  type_int  i;
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  /*Lee archivos de la simulacion*/
  read_gadget();

  for(i=0;i<3;i++)
  {
    pmin[i] = 1.E26; 
    pmax[i] = -1.E26;
  }

  for(i=0;i<cp.npart;i++)
  {
    if(P.x[i] > pmax[0]) pmax[0] = P.x[i];
    if(P.x[i] < pmin[0]) pmin[0] = P.x[i];
    if(P.y[i] > pmax[1]) pmax[1] = P.y[i];
    if(P.y[i] < pmin[1]) pmin[1] = P.y[i];
    if(P.z[i] > pmax[2]) pmax[2] = P.z[i];
    if(P.z[i] < pmin[2]) pmin[2] = P.z[i];
  }

  printf("xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  printf("ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  printf("zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);
 
  for(i=0;i<cp.npart;i++)
  {
    P.x[i] -= pmin[0];
    P.y[i] -= pmin[1];
    P.z[i] -= pmin[2];
  }

  cp.lbox = pmax[0]-pmin[0];
  for(i = 1; i < 3; i++)
    if(cp.lbox < (pmax[i] - pmin[i])) cp.lbox = (pmax[i] - pmin[i]);

  cp.lbox *= 1.001;

  //cp.lbox *= POSFACTOR;
  //cp.lbox *= 1.001;
  
  fprintf(stdout,"cp.lbox %f....\n",cp.lbox);
  GREEN("Fin Change Positions\n");
  fflush(stdout);
 
  iden.r0 = (double *) malloc(nfrac*sizeof(double));

  for(i=0;i<nfrac;i++)
  {
    fprintf(stdout, "\nBegins Identification : Step %d of %d \n",i+1,nfrac);
    
    iden.r0[i]  = fof[i];
    iden.r0[i] *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0;

    if(iden.r0[i] <= cp.soft)
    {
      fprintf(stdout,"cambia Linking length = %f \n",iden.r0[i]);
      iden.r0[i] = cp.soft;
    }

    fprintf(stdout,"Linking length = %f \n",iden.r0[i]);
  }

  iden.nobj = cp.npart;

  fprintf(stdout, "\nBegins Identification\n");
  identification();

  /************* TERMINO LA IDENTIFICACION ***************/

  free_particles();
  grid_free();

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}
