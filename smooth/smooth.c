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
#include "libfitpack.h"
#ifdef PARTICLES
  #include "sph.h"
#endif

void write(const type_int NNN, const type_real * restrict fof);

int main(int argc, char **argv)
{
  type_int    i,NNN;
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  NNN  = atoi(argv[2]);

  //Lee archivos de la simulacion //
  leeheader();

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  fprintf(stdout,"lbox %g Kpc\n",cp.lbox);
  GREEN("**********************************\n");

  read_grup_fof(fof);
  read_segment(NNN,fof);

  suaviza();

  #ifdef PARTICLES
    calculo_rho();
  #endif

  write(NNN,fof);

  for(i=0;i<cp.nseg;i++)
    free(Seg[i].Pos_list);
  free(Seg);
  free(Gr);

  TIMER(end);

  fprintf(stdout,"Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}

static void set_name(const char * prefix, char * name, const type_int NNN, const type_real * fof)
{

  sprintf(name,"%.2d_%.4d_%s",snap.num,NNN,prefix);

  #ifdef MCRITIC
	  sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);

	return;
}

void write(const type_int NNN, const type_real * restrict fof)
{

  type_int i,j;
  char filename[200];
  FILE *pfout, *pfpropiedades;

	#ifdef NEW
		set_name("smooth_new_segmentos",filename,NNN,fof);
	#else
		set_name("smooth_segmentos",filename,NNN,fof);
	#endif
  pfout=fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pfout);

	#ifdef NEW
  	set_name("smooth_new_propiedades",filename,NNN,fof);
	#else
	  set_name("smooth_propiedades",filename,NNN,fof);
	#endif
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
    #ifdef PARTICLES
      fwrite(&Seg[i].mass_part,sizeof(type_real),1,pfpropiedades);
      fwrite(&Seg[i].vol,sizeof(type_real),1,pfpropiedades);
      fwrite(&Seg[i].rho,sizeof(type_real),1,pfpropiedades);
      fwrite(&Seg[i].mu,sizeof(type_real),1,pfpropiedades);
    #endif

    fwrite(&Seg[i].size,sizeof(type_int),1,pfout);
    for(j=0;j<Seg[i].size;j++)
    {
      fwrite(&Seg[i].Pos_list[3*j+0],sizeof(type_real),1,pfout); 
      fwrite(&Seg[i].Pos_list[3*j+1],sizeof(type_real),1,pfout); 
      fwrite(&Seg[i].Pos_list[3*j+2],sizeof(type_real),1,pfout); 
    }
  }

  fclose(pfout);
  fclose(pfpropiedades);

  return;

}
