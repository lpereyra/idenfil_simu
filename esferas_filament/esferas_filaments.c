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
#include "grid.h"

int main(int argc, char **argv)
{
  double start,end;
  char filename[200];
  FILE *pfout, *pfpropiedades;
  type_int i,j,NNN,flag,itera;
  type_real r_aux = RAUX;   // TAMAÑO DEL CELL

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  #ifdef PARTICLES
  read_gadget();
  #else
  leeheader();
  #endif
  
  #if defined(ORIGINAL) || !defined(PARTICLES)
    read_grup_fof(fof);
  #endif
  read_segment(NNN,fof);

  #if defined(ORIGINAL) && defined(PARTICLES)
  free(Gr);
  #endif

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(message,"lbox %g Kpc\n",cp.lbox);GREEN(message);
  GREEN("*********************************\n");

  #ifdef PARTICLES
    grid.nobj = cp.npart;
  #else
    grid.nobj = cp.ngrup;
  #endif
  grid.ngrid = (long)(cp.lbox/r_aux);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }else{
    fprintf(stdout,"Using NGRID = %lu\n",grid.ngrid);
  }

  fflush(stdout);

  grid_init();
  grid_build();

  GREEN("******************************\n");
  BLUE("************ INIT ************\n");
  #ifdef SPH
    BLUE("************ SPH *************\n");
    sprintf(message,"******* SOFT  %g Kpc *******\n",H_SOFT);BLUE(message);
    assert(r_aux>=3.0*H_SOFT);
  #else
    BLUE("*********** K NEIGH **********\n");
    sprintf(message,"*********** NEIGH %d *********\n",K_NEIGH_SIZE);
    BLUE(message);
  #endif
  GREEN("******************************\n");
  fflush(stdout);

  itera = 1;
  flag = cp.ncent_seg;

  while(flag!=0)
  {
    #ifdef SPH
      propiedades(NNN,fof,3.0*H_SOFT,itera,flag);
    #else
      propiedades(NNN,fof,r_aux,itera,flag);
    #endif      
    itera++;
  }

  set_name(filename,NNN,"segmentos");
  pfout = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pfout);

  set_name(filename,NNN,"propiedades");
  pfpropiedades=fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pfpropiedades);

  for(i=0;i<cp.nseg;i++)
  {   
    fwrite(&Seg[i].flag,sizeof(type_int),1,pfpropiedades);
    fwrite(&Seg[i].size,sizeof(type_int),1,pfpropiedades);
    fwrite(&Seg[i].Mass[0],sizeof(type_real),2,pfpropiedades); 
    fwrite(&Seg[i].Vnodos[0],sizeof(type_real),6,pfpropiedades); 
    fwrite(&Seg[i].razon,sizeof(type_real),1,pfpropiedades); 
    fwrite(&Seg[i].len,sizeof(type_real),1,pfpropiedades);
    fwrite(&Seg[i].elong,sizeof(type_real),1,pfpropiedades);
    fwrite(&Seg[i].rms,sizeof(type_real),1,pfpropiedades);

    fwrite(&Seg[i].size,sizeof(type_int),1,pfout);

    for(j=0;j<Seg[i].size;j++)
    {
      fwrite(&Seg_centr[Seg[i].start+j].Pos[0],sizeof(type_real),1,pfout); 
      fwrite(&Seg_centr[Seg[i].start+j].Pos[1],sizeof(type_real),1,pfout); 
      fwrite(&Seg_centr[Seg[i].start+j].Pos[2],sizeof(type_real),1,pfout); 
      fwrite(&Seg[i].dens[j],sizeof(type_real),1,pfout); 
    }

  }

  fclose(pfout);
  fclose(pfpropiedades);

  grid_free();
  free(Seg);
  free(Seg_centr);
  #if defined(ORIGINAL) || !defined(PARTICLES)
  free(Gr);
  #endif
  #if defined(PARTICLES)
  free(P);
  #endif 

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}
