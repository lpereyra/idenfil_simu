#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

struct io_header header;
struct particle_data *P;
type_int *Index;
struct grup_data *Gr;
type_real pmin[3], pmax[3];

void read_gadget(void)
{
  char filename[200];
  int  ifile,ind;
  size_t total_memory;
  #ifdef MPC
  type_real POSFACTOR = 1000.;
  #endif

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
  total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  printf("Allocating %.5zu Gb for %d particles\n",total_memory,cp.npart);
  P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  Index = (type_int *) malloc((cp.npart+1)*sizeof(type_int)); // por si los id arrancan en 1
  assert(P != NULL);

  #ifdef MPC
  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(filename,"Reescala lbox %g Mpc to %g Kpc\n",cp.lbox/POSFACTOR,cp.lbox);GREEN(filename);
  GREEN("**********************************\n");
  #endif
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
  type_int d1,d2;

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

  // Definicion estructura cosmoparam
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
  printf("  Softening = %g\n",cp.soft);
  printf("*********************************** \n");
  printf("*********************************** \n");
}

void lee(char *filename, struct particle_data *Q, int *ind){
  FILE *pf;
  type_int d1, d2;
  type_int k, pc, n;
  #ifdef MPC
  type_real POSFACTOR = 1000.;
  #endif

  type_int id;
  type_real r[3];
  #ifdef STORE_VELOCITIES
  type_real v[3];
  #endif

  for(k = 0; k < 3; k++)
  {
    pmin[k] = 1.E26; 
    pmax[k] = -1.E26;
  }

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
            
        if(Q[*ind+pc].Pos[0] > pmax[0]) pmax[0] = Q[*ind+pc].Pos[0];
        if(Q[*ind+pc].Pos[0] < pmin[0]) pmin[0] = Q[*ind+pc].Pos[0];
        if(Q[*ind+pc].Pos[1] > pmax[1]) pmax[1] = Q[*ind+pc].Pos[1];
        if(Q[*ind+pc].Pos[1] < pmin[1]) pmin[1] = Q[*ind+pc].Pos[1];
        if(Q[*ind+pc].Pos[2] > pmax[2]) pmax[2] = Q[*ind+pc].Pos[2];
        if(Q[*ind+pc].Pos[2] < pmin[2]) pmin[2] = Q[*ind+pc].Pos[2];

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
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&id, size_int, 1, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        //Q[*ind+pc].id = id;
        Index[id] = *ind+pc;
        Q[*ind+pc].sub = 0;
        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  *ind += pc;
  
  fclose(pf);
}

void select_particles_fof(type_real prefix)
{

  char filename[200];
  type_int  i,j,k,id,idv,len;
  FILE *pfin;

  sprintf(snap.root,"../");
  sprintf(filename,"%s%.2d_%.2f_fof.bin",snap.root,snap.num,prefix);

  pfin = fopen(filename,"rb"); 
  fread(&len,sizeof(int),1,pfin);

  fprintf(stdout,"Archivo de Grupos %s\n",filename);
  fflush(stdout); 

  j = 0;
  for(i=0;i<len;i++)
  {  
    // tengo que comentar esta linea y el assert para las viejas versiones
    //fread(&k,sizeof(int),1,pfin); 
    fread(&id,sizeof(int),1,pfin); 
    //assert(id==k); // debe coincidir en las versiones nuevas
   
    fread(&k,sizeof(int),1,pfin);

    j+=k;
    while(k!=0)
    {
      fread(&idv,sizeof(int),1,pfin);
      P[Index[idv]].sub = id;
      k--;
    }
  }

  fprintf(stdout,"Grupos %u Part %u\n",i,j);

  fclose(pfin);
  free(Index);

  j = 0;
  for(i=0;i<cp.npart;i++)
  {
    if(P[i].sub!=0)
    {
      P[j] = P[i];
      j++;
    }
  }

  cp.npart = j;

  P = (struct particle_data *) realloc(P,cp.npart*sizeof(struct particle_data));

  fprintf(stdout,"realoc! npart %u\n",cp.npart);

}

void read_grup_fof(type_real prefix)
{
  char filename[200];
  type_int i;
  FILE *pfin;
  #ifdef MCRITIC
  type_int j;
  FILE *pfout;
  type_real cut;
  #endif

  RED("Lectura de los grupos\n");

  sprintf(filename,"%s%.2d_%.2f_centros.bin",snap.root,snap.num,prefix);
  pfin = fopen(filename,"rb"); 

  fread(&cp.ngrup,sizeof(type_int),1,pfin);
  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

  #ifdef MCRITIC

    j = 0;

    for(i=0;i<cp.ngrup;i++)
    {
      fread(&Gr[j].save,sizeof(type_int),1,pfin);
      fread(&Gr[j].id,sizeof(type_int),1,pfin);
      fread(&Gr[j].Pos[0],sizeof(type_real),1,pfin);
      fread(&Gr[j].Pos[1],sizeof(type_real),1,pfin);
      fread(&Gr[j].Pos[2],sizeof(type_real),1,pfin);
      fread(&Gr[j].NumPart,sizeof(type_int),1,pfin);

      cut = Gr[j].NumPart*cp.Mpart; // Num de particulas por la Masa de la particula [10^10 Msol/h]
      if(cut<m_critica) continue;

      j++;
    }

    if(j==cp.ngrup)
    {
      fprintf(stdout,"Grupos %u menores a mcrit %g\n",j-cp.ngrup,m_critica*1.e10);
      RED("NO DEBERIA ESCRIBIR NINGUN ARCHIVO\n");
      fclose(pfin);
      return;
    }

    cp.ngrup = j; // reasigna el num de grupos
    Gr = (struct grup_data *) realloc(Gr,cp.ngrup*sizeof(struct grup_data));
    sprintf(filename,"%.2d_%.2f_centros_cut_%.2f.bin",snap.num,prefix,m_critica);
    pfout = fopen(filename,"w");
    fwrite(&cp.ngrup,sizeof(type_int),1,pfout);

    for(i=0;i<cp.ngrup;i++)
    {
      fwrite(&Gr[i].save,sizeof(type_int),1,pfout);
      fwrite(&Gr[i].id,sizeof(type_int),1,pfout);
      fwrite(&Gr[i].Pos[0],sizeof(type_real),1,pfout);
      fwrite(&Gr[i].Pos[1],sizeof(type_real),1,pfout);
      fwrite(&Gr[i].Pos[2],sizeof(type_real),1,pfout);
      fwrite(&Gr[i].NumPart,sizeof(type_int),1,pfout);
    }

    fclose(pfout);

    fprintf(stdout,"Grupos %u mayores a mcrit %g\n",cp.ngrup,m_critica*1.e10);
      
  #else

    for(i=0;i<cp.ngrup;i++)
    {
      fread(&Gr[i].save,sizeof(type_int),1,pfin);
      fread(&Gr[i].id,sizeof(type_int),1,pfin);
      fread(&Gr[i].Pos[0],sizeof(type_real),1,pfin);
      fread(&Gr[i].Pos[1],sizeof(type_real),1,pfin);
      fread(&Gr[i].Pos[2],sizeof(type_real),1,pfin);
      fread(&Gr[i].NumPart,sizeof(type_int),1,pfin);
    }

    fprintf(stdout,"Grupos %u\n",cp.ngrup);
 
  #endif

  RED("Finaliza la lectura de los grupos\n"); fflush(stdout);

  fclose(pfin);

}
