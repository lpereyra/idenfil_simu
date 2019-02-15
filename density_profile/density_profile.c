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
  type_int    i,NNN;
  double start,end;
  type_int j;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  #ifdef ORIGINAL
  read_grup_fof(fof);
  #endif
  read_segment(NNN,fof);

  j = 0;
  for(i=0;i<cp.nseg;i++)
  {
    if(Seg[i].flag != TYPE_FLAG) continue;
    #ifdef CUT_IN_LEN
    if(Seg[i].len<LEN_MIN || Seg[i].len>LEN_MAX) continue;
    #endif
    Seg[j] = Seg[i];
    j++;
  }

  Seg = (struct segmentstd *) realloc(Seg,j*sizeof(struct segmentstd));
  cp.nseg = j;

  GREEN("********** IMPORTANTE ***********\n");
  #ifdef CUT_IN_LEN
    sprintf(message,"cut FLAG == %d && \
                     LEN MIN %f Mpc & LEN MAX %f REALOCATEA %d fil\n", \
                     TYPE_FLAG,LEN_MIN/1000.0f,LEN_MAX/1000.0f,cp.nseg);
  #else
    sprintf(message,"cut FLAG == %d REALOCATEA %d fil\n",TYPE_FLAG,cp.nseg);
  #endif
  GREEN(message);
  GREEN("**********************************\n");
  fflush(stdout);

  #ifdef CUT_ELONGACION
 
  j = 0;
  for(i=0;i<cp.nseg;i++)
  {
    if(Seg[i].flag != 2) continue;
    if(Seg[i].elong < CUT_ELONG) continue;

    Seg[j] = Seg[i];
    j++;
  }

  Seg = (struct segmentstd *) realloc(Seg,j*sizeof(struct segmentstd));
  cp.nseg = j;

  GREEN("********** IMPORTANTE ***********\n");
  sprintf(name,"cut ELONGACION %f REALOCATEA %d fil\n",CUT_ELONG,cp.nseg);GREEN(name);
  GREEN("**********************************\n");
  fflush(stdout);

  #endif

  pixelizado();

  //Lee archivos de la simulacion //
  read_gadget();

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  fprintf(stdout,"lbox %g Kpc\n",cp.lbox);
  GREEN("**********************************\n");

  propiedades(NNN,fof);
  /////////////////////////////////////////////////////////////////////////

  for(i=0;i<cp.nseg;i++)
  {
    free(Seg[i].list);
    #ifdef CALCULA_MEDIA
    free(Seg[i].Vmedia);
    #endif
  }
  free(Seg);
  free(Gr);
#ifndef RANDOM
  free(P);
#endif

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
