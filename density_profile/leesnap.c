#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

static void leeheader(char *filename){
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
  printf("  Softening = %g\n",cp.soft);
  printf("*********************************** \n");
  printf("*********************************** \n");

}

static void lee(char *filename, struct particle_data *Q, type_int *ind){
  FILE *pf;
  type_int d1, d2;
  type_int k, pc, n;
  type_real r[3];
#ifdef STORE_VELOCITIES
  type_real v[3];
#endif
#ifdef STORE_IDS
  type_int id;
#endif
 
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
      fread(&r[0], sizeof(type_real), 3, pf);
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
      fread(&v[0], sizeof(type_real), 3, pf);
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

extern void read_gadget(void)
{
  char filename[200];
  type_int  ifile,ind;

  ifile = snap.nfiles;
  if(ifile>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  size_t total_memory;

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
  total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  printf("Allocating %.5zu Gb for %u particles\n",total_memory,cp.npart);
  P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  assert(P != NULL);
  
  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++)
  {
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,P,&ind);
  }

  fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);

  return;
}

static void set_name(const char * prefix, char * name, const type_int NNN, const type_real * fof)
{
  sprintf(name,"../smooth_particles_sigma_2_rvir/");
  sprintf(name,"%s%.2d_%.4d",name,snap.num,NNN);
  sprintf(name,"%s_smooth",name);
  sprintf(name,"%s_%s",name,prefix);
  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);  

	return;
}

void read_segment(type_int NNN, type_real *fof)
{
  char  filename[200];
  type_int i,j,k,l,c;
  type_int *aux_flag;
  type_int  flag;
  type_real Mass[2];
  type_real Vnodos[6];
  type_real razon;
  type_real len;
  type_real elong;
  type_real rms;
  type_real aux_len;
  type_real mass_part;
  type_real vol;
  type_real rho;
  type_real mu;
  type_real sigma;
#ifdef CUT_ELONGACION
  type_real r[3];
#endif
  FILE  *pf;

  aux_len  = cbrt(3.0/(4.0*M_PI*cp.Mpart*(type_real)cp.npart));
  aux_len *= cp.lbox;
  aux_len *= fof[1];

  set_name("propiedades",filename,NNN,fof);
  fprintf(stdout,"%s\n",filename);

  pf = fopen(filename,"rb"); 

  fread(&cp.nseg,sizeof(type_int),1,pf);

  Seg = (struct segmentstd *) malloc(cp.nseg*sizeof(struct segmentstd));
  aux_flag = (type_int *) malloc(cp.nseg*sizeof(type_int));

  fprintf(stdout,"Propiedades Segmentos %d\n",cp.nseg);
  fflush(stdout);

  cp.ngrup = j =0;
  for(i=0;i<cp.nseg;i++)
  { 
    aux_flag[i] = 0;

    Seg[j].id = i;

    fread(&flag,sizeof(type_int),1,pf);
    fread(&Seg[j].size,sizeof(type_int),1,pf);
    fread(Mass,sizeof(type_real),2,pf);
    fread(Vnodos,sizeof(type_real),6,pf);
    fread(&razon,sizeof(type_real),1,pf);
    fread(&len,sizeof(type_real),1,pf);
    fread(&elong,sizeof(type_real),1,pf);
    fread(&rms,sizeof(type_real),1,pf);
    fread(&mass_part,sizeof(type_real),1,pf);
    fread(&vol,sizeof(type_real),1,pf);
    fread(&rho,sizeof(type_real),1,pf);
    fread(&mu,sizeof(type_real),1,pf);
    fread(&sigma,sizeof(type_real),1,pf);

    Seg[j].Rvir_2[0] = RVIR_FACTOR*cbrt(Mass[0])*aux_len;
    Seg[j].Rvir_2[1] = RVIR_FACTOR*cbrt(Mass[1])*aux_len;

    Seg[j].Rvir_2[0] *= Seg[j].Rvir_2[0];
    Seg[j].Rvir_2[1] *= Seg[j].Rvir_2[1];

    if(flag != TYPE_FLAG) continue;
    #ifdef CUT_IN_LEN
    if(len<LEN_MIN || len>LEN_MAX) continue;
    #endif

    cp.ngrup += Seg[j].size;
    aux_flag[i] = 1;
    j++;
  }

  fclose(pf);

  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));
  Seg = (struct segmentstd *) realloc(Seg,j*sizeof(struct segmentstd));
  cp.nseg = j;

  set_name("segmentos",filename,NNN,fof);

  fprintf(stdout,"%s\n",filename);

  pf = fopen(filename,"rb"); 
  fread(&k,sizeof(type_int),1,pf);

  j = c = 0;
  for(i=0;i<k;i++)
  {
    if(aux_flag[i]==0)
    {

      fread(&l,sizeof(type_int),1,pf);
      fseek(pf,3*l*sizeof(type_real),SEEK_CUR);

    }else{

      fread(&l,sizeof(type_int),1,pf);

      assert(l==Seg[j].size);

      Seg[j].start = c;

      for(l=0;l<Seg[j].size;l++)
      {
        fread(Gr[c].Pos, sizeof(type_real), 3, pf);
        c++;
      }

      j++;
    }
  }

  fclose(pf);
  free(aux_flag);

  assert(j==cp.nseg);

  GREEN("********** IMPORTANTE ***********\n");
  #ifdef CUT_IN_LEN
    sprintf(message,"cut FLAG == %d && \
                     LEN MIN %f Mpc & LEN MAX %f REALOCATEA %d fil\n", \
                     TYPE_FLAG,LEN_MIN/1000.0f,LEN_MAX/1000.0f,cp.nseg);
  #else
    sprintf(message,"cut FLAG == %d REALOCATEA %d fil\n",TYPE_FLAG,cp.nseg);
  #endif
  GREEN(message);

#ifdef CUT_ELONGACION

  l = c = 0;

  for(i=0;i<cp.nseg;i++)
  {
    len = 0.0f;
    for(k=Seg[i].start+1;k<Seg[i].start+Seg[i].size;k++)
    {
      for(j=0;j<3;j++)
      {
        r[j] = Gr[k].Pos[j]-Gr[k-1].Pos[j];    
        #ifdef PERIODIC
        if(r[j]> 0.5*cp.lbox) r[j] -= cp.lbox;
        if(r[j]<-0.5*cp.lbox) r[j] += cp.lbox;
        #endif
      }
      len += sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    }

    for(j=0;j<3;j++)
    {
      r[j] = Gr[Seg[i].start+Seg[i].size-1].Pos[j]-Gr[Seg[i].start].Pos[j];    
      #ifdef PERIODIC
      if(r[j]> 0.5*cp.lbox) r[j] -= cp.lbox;
      if(r[j]<-0.5*cp.lbox) r[j] += cp.lbox;
      #endif
    }
    
    elong = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);

    elong /= len;
    
    if(elong<CUT_ELONG)
      continue;

    for(k=Seg[i].start;k<Seg[i].start+Seg[i].size;k++)
    {
      Gr[c].Pos[0] = Gr[k].Pos[0];
      Gr[c].Pos[1] = Gr[k].Pos[1];
      Gr[c].Pos[2] = Gr[k].Pos[2];
      c++;
    }

    Seg[i].start = c-Seg[i].size;
    Seg[l] = Seg[i];

    l++;

  }

  cp.nseg = l;
  cp.ngrup = c;

  Seg = (struct segmentstd *) realloc(Seg,cp.nseg*sizeof(struct segmentstd));
  Gr  = (struct grup_data  *) realloc(Gr,cp.ngrup*sizeof(struct grup_data));

  fprintf(stdout,"Segmentos CUT ELONGATION %d\n",cp.nseg);
  fprintf(stdout,"Nodos     CUT ELONGATION %d\n",cp.ngrup);

#endif

  return;

}
