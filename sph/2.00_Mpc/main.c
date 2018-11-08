#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "timer.h"
#include "colores.h"
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"

static void Write();

int main(int argc, char **argv)
{
  type_int NNN;
  double start,end;

  TIMER(start);

  init_variables(argc,argv);
  omp_set_nested(1);
 
  // Lee archivos de la simulacion 
  NNN  = atoi(argv[2]);

  read_segment(NNN, fof);

  Write();

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}

static void Write()
{
  type_int i;
  char filename[200];
  FILE *pfout;
  #ifdef WRITE_ASCII
   FILE *pfout_ascii;
  #endif
 
  sprintf(filename,"mu_%.2f_Mpc.bin",HSPH*1e-3);
  pfout = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pfout);          
  #ifdef WRITE_ASCII
    sprintf(filename,"mu_%.2f_Mpc.dat",HSPH*1e-3);
    pfout_ascii = fopen(filename,"w");
    fprintf(pfout_ascii,"%d\n",cp.npart);         
  #endif

  for(i=0;i<cp.nseg;i++)
  {
    fwrite(&Seg[i].id,sizeof(type_int),1,pfout);          
    fwrite(&Seg[i].flag,sizeof(type_int),1,pfout);          
    fwrite(&Seg[i].rho,sizeof(type_real),1,pfout);         
    fwrite(&Seg[i].mu,sizeof(type_real),1,pfout);         
    fwrite(&Seg[i].mass_part,sizeof(type_real),1,pfout);
    fwrite(&Seg[i].vol,sizeof(type_real),1,pfout);         
    fwrite(Seg[i].Mass,sizeof(type_real),2,pfout);         
    fwrite(Seg[i].Rvir,sizeof(type_real),2,pfout);         
    fwrite(&Seg[i].razon,sizeof(type_real),1,pfout);         
    fwrite(&Seg[i].len,sizeof(type_real),1,pfout);         
    fwrite(&Seg[i].elong,sizeof(type_real),1,pfout);         

    #ifdef WRITE_ASCII
      fprintf(pfout_ascii,"%d %f %f\n",Seg[i].id,Seg[i].rho,Seg[i].mass_part);
    #endif
  }

  fclose(pfout);
  #ifdef WRITE_ASCII
   fclose(pfout_ascii);
  #endif

  return;
}
