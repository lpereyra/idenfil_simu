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

static void cut_fil(void);

int main(int argc, char **argv)
{
  type_int i, NNN;
  double start,end;
  char filename[200];
  FILE *pf;

  TIMER(start);
  
  init_variables(argc,argv);

  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  leeheader();
  
  read_grup_fof(fof);
  read_segment(NNN,fof);
  cut_fil();

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(message,"lbox %g Kpc\n",cp.lbox);GREEN(message);
  GREEN("*********************************\n");

  #ifdef CENTERS_RANDOM
    sprintf(message,"CREA RANDOMS DE ACUERDO A LA MASA DE LAS PUNTAS\n");RED(message);
    crea_random(NNN);
    sprintf(message,"FINALIZA\n");GREEN(message);
  #endif

  #ifdef CENTERS_RANDOM
    sprintf(filename,"%.2d_%.4d_pares_random_halos_%.2f_%.2f.bin",snap.num,NNN,fof[0],fof[1]);
  #else
    sprintf(filename,"%.2d_%.4d_pares_halos_%.2f_%.2f.bin",snap.num,NNN,fof[0],fof[1]);
  #endif

  pf = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(int),1,pf);
  fprintf(stdout,"Total de Pares %d\n",cp.nseg); 

  for(i=0;i<cp.nseg;i++)
  {  
    fwrite(&Gr[Seg[i].list[0]].id,sizeof(type_int),1,pf);
    fwrite(&Gr[Seg[i].list[1]].id,sizeof(type_int),1,pf);
    fwrite(&Seg[i].len,sizeof(type_real),1,pf);
    fwrite(&Gr[Seg[i].list[0]].Mass,sizeof(type_real),1,pf);
    fwrite(&Gr[Seg[i].list[1]].Mass,sizeof(type_real),1,pf);
    fwrite(&Gr[Seg[i].list[0]].vcm,sizeof(type_real),3,pf);
    fwrite(&Gr[Seg[i].list[1]].vcm,sizeof(type_real),3,pf);
  }

  fclose(pf);

  free(Seg);
  free(Gr);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}

static void cut_fil(void)
{
  type_int i, j;

  j = 0;
  for(i=0;i<cp.nseg;i++)
  {
    if(Seg[i].flag != 2) continue;
    #ifdef CUT_IN_LEN
    if(Seg[i].len<LEN_MIN || Seg[i].len>LEN_MAX) continue;
    #endif

    Seg[j] = Seg[i];
    j++;
  }

  Seg = (struct segmentstd *) realloc(Seg,j*sizeof(struct segmentstd));
  cp.nseg = j;

  #ifdef CUT_IN_LEN
  GREEN("********** IMPORTANTE ***********\n");
  sprintf(message,"cut LEN %f Mpc REALOCATEA %d fil\n",0.5*(LEN_MIN+LEN_MAX)/1000.0f,cp.nseg);GREEN(message);
  GREEN("**********************************\n");
  fflush(stdout);
  #else
  GREEN("********** IMPORTANTE ***********\n");
  sprintf(message,"only TYPE 2 Mpc REALOCATEA %d fil\n",cp.nseg);GREEN(message);
  GREEN("**********************************\n");
  fflush(stdout);
  #endif


  return;
}

