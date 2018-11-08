#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
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

  type_real r[3];
  #ifdef STORE_VELOCITIES
  type_real v[3];
  #endif
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
        Q[*ind+pc].Pos[0] = r[0]*POSFACTOR;
        Q[*ind+pc].Pos[1] = r[1]*POSFACTOR;
        Q[*ind+pc].Pos[2] = r[2]*POSFACTOR;
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

void read_grup_fof(type_real prefix)
{
  char filename[200];
  int  i;
  int j;
  FILE *pfin;
  #ifdef MCRITIC
  FILE *pfout;
  type_real cut;
  #endif

  RED("Lectura de los grupos\n");

  sprintf(snap.root,"./");
  sprintf(filename,"%s%.2d_%.2f_centros.bin",snap.root,snap.num,prefix);
  pfin = fopen(filename,"rb"); 

  fread(&cp.ngrup,sizeof(int),1,pfin);
  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));
     
  j = 0;
  for(i=0;i<cp.ngrup;i++)
  {
    fread(&Gr[j].save,sizeof(int),1,pfin);
    fread(&Gr[j].id,sizeof(type_int),1,pfin);
    fread(&Gr[j].Pos[0],sizeof(type_real),1,pfin);
    fread(&Gr[j].Pos[1],sizeof(type_real),1,pfin);
    fread(&Gr[j].Pos[2],sizeof(type_real),1,pfin);
    fread(&Gr[j].NumPart,sizeof(int),1,pfin);
    #ifdef MCRITIC
      cut = Gr[j].NumPart*cp.Mpart; // Num de particulas por la Masa de la particula [10^10 Msol/h]
      if(cut<m_critica) continue;
    #endif
    j++;
  }

  #ifdef MCRITIC

    cp.ngrup = j; // reasigna el num de grupos
    Gr = (struct grup_data *) realloc(Gr,cp.ngrup*sizeof(struct grup_data));
    fprintf(stdout,"Grupos %d mayores a mcrit %g\n",cp.ngrup,m_critica*1.e10);

  #else

    fprintf(stdout,"Grupos %d\n",cp.ngrup);

  #endif
  RED("Finaliza la lectura de los grupos\n"); fflush(stdout);

  fclose(pfin);

}

void select_particles_fof(type_real *fof)
{

  char filename[200];
  int  i,j,k,id,idv,len,tipo;
  type_int *index;
  FILE *pfin;

  index = (type_int *) malloc(cp.npart*sizeof(type_int));

  #pragma omp parallel for \
  schedule(static) default(none) private(i) shared(P,index,cp)
  for(i=0;i<cp.npart;i++)
  {
    index[P[i].id] = i;
    P[i].nfof[0] = 0;
    P[i].nfof[1] = 0;
  }

  for(tipo=0;tipo<nfrac;tipo++)
  {

    sprintf(filename,"%s%.2d_%.2f_fof.bin",snap.root,snap.num,fof[tipo]);

    pfin = fopen(filename,"rb"); 
    fread(&len,sizeof(int),1,pfin);

    fprintf(stdout,"Archivo de Grupos %s\n",filename);
    fflush(stdout); 

    j = 0;
    for(i=0;i<len;i++)
    {   
      fread(&k,sizeof(int),1,pfin);
      fread(&id,sizeof(int),1,pfin); 

      fread(&k,sizeof(int),1,pfin);

      j+=k;
      while(k!=0)
      {
        fread(&idv,sizeof(int),1,pfin);
        P[index[idv]].nfof[tipo] = id;
        k--;
      }
    }

    fprintf(stdout,"Grupos %d Part %d\n",i,j);

    fclose(pfin);

  }


  free(index);
}


