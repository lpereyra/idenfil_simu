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

void read_pares(int NNN, type_real *fof)
{
  char  filename[200];
  int   i,k;
  FILE  *pf;
 
  #ifdef MCRITIC
    sprintf(filename,"../pares_random/%.2d_%.4d_pares_halos_%.2f_%.2f_%.2f.bin",snap.num,NNN,m_critica,fof[0],fof[1]);
  #else
    sprintf(filename,"../pares_random/%.2d_%.4d_pares_halos_%.2f_%.2f.bin",snap.num,NNN,fof[0],fof[1]);
  #endif

  pf = fopen(filename,"rb"); 

  fread(&cp.nseg,sizeof(int),1,pf);

  fprintf(stdout,"Segmentos %d\n",cp.nseg);
  fflush(stdout);

  fprintf(stdout,"----------- Le resta uno cdo LEE cuidado aca\n");
  fflush(stdout);

  Seg = (struct segmentstd *) malloc(cp.nseg*sizeof(struct segmentstd));

  for(i=0;i<cp.nseg;i++)
  {
    Seg[i].size = 2;
    Seg[i].list = (int *) malloc(Seg[i].size*sizeof(int));

    for(k=0;k<Seg[i].size;k++)
    {
      fread(&Seg[i].list[k],sizeof(int),1,pf);
      Seg[i].list[k] -= 1;
      assert(Seg[i].list[k]>=0);
    }

    fread(&Seg[i].len,sizeof(type_real),1,pf);
    fread(&Seg[i].Mass[0],sizeof(type_real),1,pf);
    fread(&Seg[i].Mass[1],sizeof(type_real),1,pf);
    Seg[i].razon = Seg[i].Mass[0]/Seg[i].Mass[1];
    Seg[i].flag = 0;
    Seg[i].flag = Seg[i].Mass[0]>=NNN*cp.Mpart ? Seg[i].flag+1: Seg[i].flag;
    Seg[i].flag = Seg[i].Mass[1]>=NNN*cp.Mpart ? Seg[i].flag+1: Seg[i].flag;
  
    assert(Seg[i].flag==2);
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

#ifdef VEL_RELATIVA

void grupos_fof(type_real *fof)
{

  char filename[200];
  int  i,j,k,id,idv;
  int  size,len;
  type_int *index;
  FILE *pfin;

  index = (type_int *) malloc((cp.npart+1)*sizeof(type_int)); // por si los id arrancan en 1

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(static) default(none) private(i) shared(P,index,cp)
  for(i=0;i<cp.npart;i++)
    index[P[i].id] = i;

  sprintf(snap.root,"../../");
  sprintf(filename,"%s%.2d_%.2f_fof.bin",snap.root,snap.num,fof[1]);

  pfin = fopen(filename,"rb"); 
  fread(&len,sizeof(int),1,pfin);

  fprintf(stdout,"Archivo de Grupos %s\n",filename);
  fflush(stdout); 

  j = 0;
  for(i=0;i<len;i++)
  {   
    fread(&k,sizeof(int),1,pfin);
    fread(&id,sizeof(int),1,pfin); 
    fread(&size,sizeof(int),1,pfin);

    assert(Gr[i].NumPart==size);

    Gr[i].Vnod[0] = Gr[i].Vnod[1] = Gr[i].Vnod[2] = 0.f;

    j+=size;
    for(k=0;k<size;k++)
    {
      fread(&idv,sizeof(int),1,pfin);

      Gr[i].Vnod[0] += P[index[idv]].Vel[0];
      Gr[i].Vnod[1] += P[index[idv]].Vel[1];
      Gr[i].Vnod[2] += P[index[idv]].Vel[2];
    }

    for(k=0;k<3;k++)
      Gr[i].Vnod[k] *= (1.0f/(type_real)size);

  }

  fprintf(stdout,"Grupos %d Part %d\n",i,j);

  fclose(pfin);
  free(index);
}

#endif
