#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

void leeheader(char *filename){
  FILE *pf;
  int d1,d2;

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fclose(pf);

  /* Definicion estructura cosmoparam */
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;
  cp.npart     = header.npartTotal[1];
  cp.Mpart     = header.mass[1];
  //cp.Mpart    = 3.143E-4*cp.hparam;  /*A-5*/
  //cp.Mpart    = 3.929E-5*cp.hparam;  /*A-4*/
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

  printf("*********************************** \n");
  printf("*   Parametros de la simulacion   * \n");
  printf("*********************************** \n");
  printf("  Numero de particulas = %d \n", cp.npart);
  printf("  Lado del box = %g \n", cp.lbox);
  printf("  Redshift = %g \n", cp.redshift);
  printf("  Omega Materia = %g \n", cp.omegam);
  printf("  Omega Lambda = %g \n", cp.omegal);
  printf("  Parametro de Hubble = %g \n",cp.hparam);
  printf("  Masa por particula = %g \n",cp.Mpart);
  printf("  Softening = %g\n",cp.soft);
  printf("*********************************** \n");
  printf("*********************************** \n");
}

void read_segment(int NNN, type_real *fof)
{
  char  filename[200];
  int   i,k;
  FILE  *pf;
 
  #ifdef MCRITIC
    sprintf(filename,"../%.2d_%.4d_segmentos_cut_%.2f_%.2f_%.2f.bin",snap.num,NNN,m_critica,fof[0],fof[1]);
  #else
    sprintf(filename,"../%.2d_%.4d_segmentos_%.2f_%.2f.bin",snap.num,NNN,fof[0],fof[1]);
  #endif

  pf = fopen(filename,"rb"); 

  fread(&cp.nseg,sizeof(int),1,pf);

  fprintf(stdout,"Segmentos %d\n",cp.nseg);
  fflush(stdout);

  Seg = (struct segmentstd *) malloc(cp.nseg*sizeof(struct segmentstd));

  for(i=0;i<cp.nseg;i++)
  {
    fread(&Seg[i].size,sizeof(int),1,pf);

    Seg[i].list = (int *) malloc(Seg[i].size*sizeof(int));

    for(k=0;k<Seg[i].size;k++)
    {
      fread(&Seg[i].list[k],sizeof(int),1,pf);
    }
  }

  fclose(pf);

  #ifdef MCRITIC
  sprintf(filename,"../%.2d_%.4d_propiedades_cut_%.2f_%.2f_%.2f.bin",snap.num,NNN,m_critica,fof[0],fof[1]);
  #else
  sprintf(filename,"../%.2d_%.4d_propiedades_%.2f_%.2f.bin",snap.num,NNN,fof[0],fof[1]);
  #endif

  pf = fopen(filename,"rb"); 

  fread(&k,sizeof(int),1,pf);

  assert(k==cp.nseg);

  fprintf(stdout,"Propiedades Segmentos %d\n",cp.nseg);
  fflush(stdout);

  for(i=0;i<cp.nseg;i++)
  {  
    fread(&Seg[i].flag,sizeof(int),1,pf);
    fread(&k,sizeof(int),1,pf);
    fread(&Seg[i].Mass[0],sizeof(float),1,pf);
    fread(&Seg[i].Mass[1],sizeof(float),1,pf);
    fread(&Seg[i].razon,sizeof(float),1,pf);
    fread(&Seg[i].len,sizeof(float),1,pf);
    fread(&Seg[i].elong,sizeof(float),1,pf);
    fread(&Seg[i].rms,sizeof(float),1,pf);

    assert(k==Seg[i].size);
  }

  fclose(pf);

  return;
}

void read_grup_fof(type_real *fof)
{
  char  filename[200];
  int   i;
  FILE  *pfin;
 
  #ifdef MCRITIC
    sprintf(filename,"../%.2d_%.2f_centros_cut_%.2f.bin",snap.num,fof[1],m_critica);
  #else
    sprintf(filename,"../../%.2d_%.2f_centros.bin",snap.num,fof[1]);
  #endif

  pfin = fopen(filename,"rb"); 

  fread(&cp.ngrup,sizeof(int),1,pfin);

  fprintf(stdout,"Grupos %d\n",cp.ngrup);
  fflush(stdout);

  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

  for(i=0;i<cp.ngrup;i++)
  {
    fread(&Gr[i].save,sizeof(int),1,pfin);
    fread(&Gr[i].id,sizeof(int),1,pfin);
    fread(&Gr[i].Pos[0],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[1],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[2],sizeof(float),1,pfin);
    fread(&Gr[i].NumPart,sizeof(int),1,pfin);
  }

  fclose(pfin);

  return;

}
