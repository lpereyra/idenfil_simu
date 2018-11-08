//compilar cbuild.csh pos; ./pos

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "s2plot.h"
#include <omp.h>

#define PERIODIC   
#define POSFACTOR 1000.0
#define STORE_IDS

#include "timer.h"
#include "colores.h"
#include "cosmoparam.h"
#include "variables.c"
#include "leesnap.c"

int main(int argc, char **argv)
{
  char filename[200];
  int i,j,k;
  char xopt[] = "BCDET";
  char yopt[] = "BCDET";
  char zopt[] = "BCDET";
  double start,end;

  TIMER(start);

  init_variables(argc,argv);

  /* Lee archivos de la simulacion */
  read_gadget();
  read_grup_fof(fof[1]);
  select_particles_fof(fof);

  for(j=0;j<3;j++)
  {
    pmin[j] = 2.0f*POSFACTOR*cp.lbox;
    pmax[j] = 0.0f;
  }
 
  for(i=0;i<cp.npart;i++)
  {
    for(j=0;j<3;j++)
    {
      pmin[j] = (P[i].Pos[j]<pmin[j]) ? P[i].Pos[j] : pmin[j];
      pmax[j] = (P[i].Pos[j]>pmax[j]) ? P[i].Pos[j] : pmax[j];
    }
  }

  s2opend("/?",argc,argv);                                      // Open the display
  ss2spt(1);                                                    // Generate new projection type
  ss2sbc(0., 0., 0.);		            	                // Set background colour
  s2swin(pmin[0],pmax[0],pmin[1],pmax[1],pmin[2],pmax[2]);      // Set the window coordinates
  s2box(xopt,0,0,yopt,0,0,zopt,0,0);   		                // Draw coordinate box

  /*
  s2sci(1);				// Set the colour
  s2slw(1.0);

  srand48(15);
  for(i=0;i<cp.npart;i++)
  {
    if(drand48()>0.1) continue;
    s2pt1(P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],1);
  }
  */

  s2slw(1.0);

  for(i=0;i<cp.npart;i++)
  {
    if(P[i].nfof[0]==0) continue;
    //fprintf(stdout,"%d\n",P[i].nfof[0]);
    if(P[i].nfof[0]==1)
    {
      s2sci(4);				// Set the colour
      //s2pt1(P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],1);

      if(P[i].nfof[1]!=0)
      {
        s2sci(2);				// Set the colour
        s2pt1(P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],1);
      }
    }
  }

  s2sci(3);				// Set the colour
  s2slw(5.0);
 
  for(i=0;i<cp.ngrup;i++)
  {      
    if(Gr[i].save==1)
      s2pt1(Gr[i].Pos[0],Gr[i].Pos[1],Gr[i].Pos[2],1);
  }



  TIMER(end);
  printf("Total time %f\n",end-start);
  s2show(1);
  free(P);
  return(EXIT_SUCCESS);

}
