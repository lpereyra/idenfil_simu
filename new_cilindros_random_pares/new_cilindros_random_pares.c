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
#include "calcula.h"

static void pixelizado(void);

int main(int argc, char **argv)
{
  type_int    NNN;
  double start,end;
  #ifdef CALCULA_MEDIA
  type_int    i;
  #endif

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  read_pares(NNN,fof);
  read_grup_fof(fof);

  pixelizado();

  // Lee archivos de la simulacion //
  read_gadget();

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  fprintf(stdout,"lbox %g Kpc\n",cp.lbox);
  GREEN("**********************************\n");

  //////////////////////////////////////////////
  propiedades(NNN,fof);
  //////////////////////////////////////////////

  free(Seg_cent);
  free(Gr);
  free(P);
  #ifdef CALCULA_MEDIA
  for(i=0;i<cp.nseg;i++)
    free(Seg[i].Vmedia);
  #endif
  free(Seg);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}

static void pixelizado(void)
{
  type_int i,j,k;  
  struct grup_data *Gr_aux;
  
  j = 0;
  for(i=0;i<cp.nseg;i++)
    for(k=0;k<Seg[i].size;k++)
      j++;   

  cp.ngrup = j;
  Gr_aux = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

  j = 0;
  for(i=0;i<cp.nseg;i++)
  {
    for(k=0;k<Seg[i].size;k++)
    {      
      Gr_aux[j] = Gr[Seg[i].list[k]];
      Gr_aux[j].save = i;
      Seg[i].list[k]  = j;
      j++;   
    }
  }

  assert(j==cp.ngrup);
  sprintf(message,"Ngrup for grid %u\n",cp.ngrup);RED(message);

  free(Gr);

  //Swap arrays
  Gr     = &(* Gr_aux);
  Gr_aux = NULL;
  free(Gr_aux);

  return;

}
