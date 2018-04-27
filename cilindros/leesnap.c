#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

void read_gadget(void){
  char filename[200];
  int  ifile,ind;
  size_t total_memory;

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
  total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  printf("Allocating %.5zu Gb for %d particles\n",total_memory,cp.npart);
  P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  assert(P != NULL);

  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++){
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,P,&ind);
  }

  fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);
}

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

void lee(char *filename, struct particle_data *Q, int *ind){
  FILE *pf;
  int d1, d2;
  int k, pc, n;
  #ifdef MPC 
  type_real POSFACTOR = 1000.;
  #endif

  type_real r[3],v[3];
  type_int id;
 
  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&r[0], size_real, 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        #ifdef MPC
          Q[*ind+pc].Pos[0] = r[0]*POSFACTOR;
          Q[*ind+pc].Pos[1] = r[1]*POSFACTOR;
          Q[*ind+pc].Pos[2] = r[2]*POSFACTOR;
        #else
          Q[*ind+pc].Pos[0] = r[0];
          Q[*ind+pc].Pos[1] = r[1];
          Q[*ind+pc].Pos[2] = r[2];
        #endif
        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_VELOCITIES
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&v[0], size_real, 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        Q[*ind+pc].Vel[0] = v[0]*VELFACTOR;
        Q[*ind+pc].Vel[1] = v[1]*VELFACTOR;
        Q[*ind+pc].Vel[2] = v[2]*VELFACTOR;
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_IDS
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&id, size_int, 1, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        Q[*ind+pc].id = id;
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  *ind += pc;
  
  fclose(pf);
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
    fread(&Seg[i].razon,sizeof(float),1,pf);
    fread(&Seg[i].len,sizeof(float),1,pf);
    fread(&Seg[i].elong,sizeof(float),1,pf);
    fread(&Seg[i].rms,sizeof(float),1,pf);
    
    assert(k==Seg[i].size);
  }

  fclose(pf);

  #ifndef CALCULA_MEDIA
  
  int j,n,id,contador;

  contador = 0;
  for(i=0;i<NTHREADS;i++)
  {
    #ifdef MCRITIC
    sprintf(filename,"../%.2d_%.4d_vmedia_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NNN,m_critica,fof[0],fof[1],i);
    #else
    sprintf(filename,"../%.2d_%.4d_vmedia_%.2f_%.2f.%.2d.bin",snap.num,NNN,fof[0],fof[1],i);
    #endif

    pf = fopen(filename,"rb"); 

    fread(&k,sizeof(int),1,pf);
    fread(&n,sizeof(int),1,pf);

    for(j=0;j<k;j++)
    {
      fread(&id,sizeof(int),1,pf);
      Seg[id].Vmedia = (type_real *) malloc(3*n*sizeof(type_real));
      fread(&Seg[id].Vmedia[0],sizeof(type_real),3*n,pf);
    }
    
    contador += k;
   
    fprintf(stdout,"Tid %d Nfil %d\n",i,k);  

    fclose(pf);
  }
    
  assert(contador==cp.nseg);

  #endif

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
