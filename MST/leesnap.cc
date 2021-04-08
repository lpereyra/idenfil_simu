#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.hh"
#include "colores.hh"
#include "leesnap.hh"

static type_int *Index;
static struct io_header header;
struct particle_data *P;
struct grup_data *Gr;

static void leeheader(const char *filename)
{
  FILE *pf;
  type_int d1,d2;
#ifdef TYPE_TWO_GADGET
  type_int blocksize;
  char label[4];
#endif

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
#endif

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  assert(d1==256);

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

#ifdef TYPE_TWO_GADGET

  type_int k, n;
  type_real mass;
  int flag = 0;
  
  while(1)
  {
    fread(&d1, sizeof(d1), 1, pf);
    fread(&label, sizeof(char), 4, pf);
    fread(&blocksize, sizeof(int), 1, pf);
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);
  
    flag = strncmp(label, "MASS", 4);
    if(flag == 0)
    {
  
      fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
      fflush(stdout);
  
      fread(&d1, sizeof(d1), 1, pf);
      for(k = 0; k < N_part_types; k++){
        for(n = 0; n < header.npart[k]; n++){
          fread(&mass, sizeof(mass), 1, pf);
          if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
            cp.Mpart = (double)mass;
          }
        }
      }
      fread(&d2, sizeof(d2), 1, pf);
      assert(d1==d2);
      break;
  
    }else{
  
      fread(&d1, sizeof(d1), 1, pf);
      fseek(pf, d1,SEEK_CUR);
      fread(&d2, sizeof(d2), 1, pf);
      assert(d1==d2);
  
    }
  
  }

  sprintf(message,"change Masa por particula = %f\n",cp.Mpart);
  RED(message);

#endif

  fclose(pf);

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

	return;
}

static void lee(char *filename, type_int *ind)
{
  FILE *pf;
  type_int d1, d2;
  type_int k, pc, n;
#ifdef TYPE_TWO_GADGET
  type_int blocksize;
  char label[4];
#endif

  type_int id;
  type_real r[3];
#ifdef STORE_VELOCITIES
	type_real v[3];
#endif

	pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif
 
  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  assert(d1==256);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&r[0], sizeof(type_real), 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        P[pc].Pos[0] = r[0]*POSFACTOR;
        P[pc].Pos[1] = r[1]*POSFACTOR;
        P[pc].Pos[2] = r[2]*POSFACTOR;
        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_VELOCITIES
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&v[0], sizeof(type_real), 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        P[pc].Vel[0] = v[0]*VELFACTOR;
        P[pc].Vel[1] = v[1]*VELFACTOR;
        P[pc].Vel[2] = v[2]*VELFACTOR;
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
#endif

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&id, sizeof(type_int), 1, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
#ifdef STORE_IDS
        P[pc].Id = id;
#endif
        Index[id] = pc;
        P[pc].sub = 0;
        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  *ind = pc;
  
  fclose(pf);
}

extern void read_gadget(void)
{
  char filename[200];
	type_int ifile, ind;
  size_t total_memory;

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
  total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  printf("Allocating %.5zu Gb for %u particles\n",total_memory,cp.npart);
  P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  assert(P != NULL);
  Index = (type_int *) calloc((cp.npart+1),sizeof(type_int)); // por si los id arrancan en 1
  assert(Index != NULL);
  
  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(filename,"Reescala\n");GREEN(filename);
  GREEN("**********************************\n");
  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++){
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,&ind);
  }

  fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);
}

extern void select_particles_fof(type_int idx)
{

  char filename[200];
  type_int  i,j,k,id,idv,len;
  FILE *pfin;
	type_real prefix = fof[idx];

  sprintf(snap.root,"../mendieta/");
  sprintf(filename,"%s%.2d_%.2d_level_%.2f_fof.bin",snap.root,snap.num,idx,prefix);

  fprintf(stdout,"Archivo de Grupos %s\n",filename);
  fflush(stdout); 

  pfin = fopen(filename,"rb"); 
  fread(&len,sizeof(int),1,pfin);

	fprintf(stdout,"ngrups %d\n",len);
  fflush(stdout); 

  j = 0;
  for(i=0;i<len;i++)
  {  
    // tengo que comentar esta linea y el assert para las viejas versiones
    fread(&k,sizeof(type_int),1,pfin); 
    fread(&id,sizeof(type_int),1,pfin); 
   
    fread(&k,sizeof(type_int),1,pfin);

    j+=k;
    while(k!=0)
    {
      fread(&idv,sizeof(type_int),1,pfin);
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

	return;
}

extern void read_grup_fof(type_real prefix)
{
  char filename[200];
  type_int i,jj;
  FILE *pfin;

  RED("Lectura de los grupos\n");

  sprintf(filename,"%s%.2d_01_level_%.2f_centros.bin",snap.root,snap.num,prefix);
  pfin = fopen(filename,"rb"); 

  fread(&cp.ngrup,sizeof(type_int),1,pfin);
  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

  for(i=0;i<cp.ngrup;i++)
  {
    fread(&Gr[i].save,sizeof(type_int),1,pfin);
    fread(&Gr[i].id,sizeof(type_int),1,pfin);
    fread(&Gr[i].NumPart,sizeof(type_int),1,pfin);
    fread(&Gr[i].Pos[0],sizeof(type_real),3,pfin);
#ifdef STORE_VELOCITIES  
	  fread(&Gr[i].vcm[0],sizeof(type_real),3,pfin);
#ifdef COMPUTE_EP
	  fread(&Gr[i].mostbound[0],sizeof(type_real),3,pfin);
#endif
	  fread(&Gr[i].sig[0],sizeof(type_real),3,pfin);
	  fread(&Gr[i].L[0],sizeof(type_real),3,pfin);
#ifdef COMPUTE_EP
	  fread(&Gr[i].lambda,sizeof(type_real),1,pfin);
#endif
	  fread(&Gr[i].m200,sizeof(type_real),1,pfin);
	  fread(&Gr[i].r200,sizeof(type_real),1,pfin);
	  fread(&Gr[i].v200,sizeof(type_real),1,pfin);
	  fread(&Gr[i].mvir,sizeof(type_real),1,pfin);
	  fread(&Gr[i].rvir,sizeof(type_real),1,pfin);
	  fread(&Gr[i].vvir,sizeof(type_real),1,pfin);
	  fread(&Gr[i].vmax,sizeof(type_real),1,pfin);
#ifdef COMPUTE_EP
	  fread(&Gr[i].Ep,sizeof(type_real),1,pfin);
	  fread(&Gr[i].Ec,sizeof(type_real),1,pfin);
#endif  
#endif  
	  fread(&Gr[i].aa,sizeof(type_real),1,pfin);
	  fread(&Gr[i].bb,sizeof(type_real),1,pfin);
	  fread(&Gr[i].cc,sizeof(type_real),1,pfin);
	  for(jj=0;jj<3;jj++)
	    fread(&Gr[i].evec[jj][0],sizeof(type_real),3,pfin);
#ifdef STORE_VELOCITIES  
	  fread(&Gr[i].aa_vel,sizeof(type_real),1,pfin);
	  fread(&Gr[i].bb_vel,sizeof(type_real),1,pfin);
	  fread(&Gr[i].cc_vel,sizeof(type_real),1,pfin);
	  for(jj=0;jj<3;jj++)
	    fread(&Gr[i].evec_vel[jj][0],sizeof(type_real),3,pfin);
#endif  
    }

    fprintf(stdout,"Grupos %u\n",cp.ngrup);
 
  RED("Finaliza la lectura de los grupos\n"); fflush(stdout);

  fclose(pfin);

	return;

}
