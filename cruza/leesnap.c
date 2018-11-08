#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

extern void leeheader(void){
  FILE *pf;
  type_int d1,d2;
  char  filename[200];

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

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
  printf("  Numero de particulas = %u \n", cp.npart);
  printf("  Lado del box = %g \n", cp.lbox);
  printf("  Redshift = %g \n", cp.redshift);
  printf("  Omega Materia = %g \n", cp.omegam);
  printf("  Omega Lambda = %g \n", cp.omegal);
  printf("  Parametro de Hubble = %g \n",cp.hparam);
  printf("  Masa por particula = %g \n",cp.Mpart);
  printf("*********************************** \n");
  printf("*********************************** \n");
}

extern void read_segment(type_int flag_cut)
{
  char  filename[200], chain[200], *p;
  struct segmentstd *Aux_Seg;
  type_int   id,  i, k;
  FILE  *pf;
 
  if(flag_cut == 0)
    sprintf(filename,"%s%s",fil.froot,fil.fname);
  else
    sprintf(filename,"%s%s",fil_cut.froot,fil_cut.fname);

  fprintf(stdout,"%s\n",filename);

  pf = fopen(filename,"rb"); 

  fread(&i,sizeof(type_int),1,pf);

  fprintf(stdout,"Segmentos %d\n",i);
  fflush(stdout);

  if(flag_cut == 0)
    cp.nseg     = i;
  else
    cp.nseg_cut = i;

  Aux_Seg = (struct segmentstd *) malloc(i*sizeof(struct segmentstd));

  for(id=0;id<i;id++)
  {
    fread(&Aux_Seg[id].size,sizeof(type_int),1,pf);
    Aux_Seg[id].list = (type_int *) malloc(Aux_Seg[id].size*sizeof(type_int));

    for(k=0;k<Aux_Seg[id].size;k++)
    {
      fread(&(Aux_Seg[id].list[k]),sizeof(type_int),1,pf);
    }
  }

  fclose(pf);

  ////////////////////////////////////////////////////////////////////

  p = strstr(filename, "segmentos");

  // Copy the part of the old filename *before* the replacement
  memcpy(chain, filename, p - filename);

  // Copy the replacement
  memcpy(chain + (p - filename), "propiedades", strlen("propiedades"));

  // Copy the rest
  strcpy(chain + (p - filename) + strlen("propiedades"), p + strlen("segmentos"));

  fprintf(stdout,"%s\n",chain);

  ////////////////////////////////////////////////////////////////////

  pf = fopen(chain,"rb"); 

  fread(&k,sizeof(type_int),1,pf);

  assert(k==i);

  fprintf(stdout,"Propiedades Segmentos %d\n",i);
  fflush(stdout);

  for(id=0;id<i;id++)
  {  
    fread(&Aux_Seg[id].flag,sizeof(type_int),1,pf);
    fread(&k,sizeof(type_int),1,pf);
    fread(&Aux_Seg[id].Mass[0],sizeof(float),2,pf);
    fread(&Aux_Seg[id].Vnodos[0],sizeof(float),6,pf);
    fread(&Aux_Seg[id].razon,sizeof(float),1,pf);
    fread(&Aux_Seg[id].len,sizeof(float),1,pf);
    fread(&Aux_Seg[id].elong,sizeof(float),1,pf);
    fread(&Aux_Seg[id].rms,sizeof(float),1,pf);
    
    assert(k==Aux_Seg[id].size);
  }

  fclose(pf);

  if(flag_cut == 0)
    Seg = &(* Aux_Seg);
  else
    Seg_cut = &(* Aux_Seg);

  Aux_Seg = NULL;
  free(Aux_Seg);

  return;
}

extern void read_grup_fof(type_int flag_cut)
{
  char  filename[200];
  struct grup_data *Aux_Gr;
  type_int   i,id;
  FILE  *pfin;

  if(flag_cut == 0)
    sprintf(filename,"%s%s",fil.croot,fil.cname);
  else
    sprintf(filename,"%s%s",fil_cut.croot,fil_cut.cname);

  pfin = fopen(filename,"rb"); 

  fread(&i,sizeof(type_int),1,pfin);

  if(flag_cut == 0)
    cp.ngrup     = i;
  else
    cp.ngrup_cut = i;

  fprintf(stdout,"Grupos %d\n",i);
  fflush(stdout);

  Aux_Gr = (struct grup_data *) malloc(i*sizeof(struct grup_data));

  for(id=0;id<i;id++)
  {
    fread(&Aux_Gr[id].save,sizeof(type_int),1,pfin);
    fread(&Aux_Gr[id].id,sizeof(type_int),1,pfin);
    fread(&Aux_Gr[id].Pos[0],sizeof(float),1,pfin);
    fread(&Aux_Gr[id].Pos[1],sizeof(float),1,pfin);
    fread(&Aux_Gr[id].Pos[2],sizeof(float),1,pfin);
    fread(&Aux_Gr[id].NumPart,sizeof(type_int),1,pfin);
  }

  fclose(pfin);

  if(flag_cut == 0)
    Gr = &(* Aux_Gr);
  else
    Gr_cut = &(* Aux_Gr);

  Aux_Gr = NULL;
  free(Aux_Gr);

  return;

}
