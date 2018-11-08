#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"
#include "grid.h"
#include "sph.h"

static struct io_header header;  

static void leeheader(char *filename)
{
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

  // Definicion estructura cosmoparam //
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;
  cp.npart     = header.npartTotal[1];
  cp.Mpart     = header.mass[1];
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

  /*
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
  printf("*********************************** \n");
  printf("*********************************** \n");
  */
}

static void lee(char *filename, struct particle_data *Q, type_int *ind)
{
  FILE *pf;
  int d1, d2;
  int k, pc, n;

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

  //fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = 0; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&r[0], sizeof(type_int), 3, pf);
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
      fread(&id, sizeof(type_int), 1, pf);
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

static void read_gadget(char *filename)
{
  type_int ind;
  //size_t total_memory;

  leeheader(filename);

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/

  cp.lbox *= POSFACTOR;
  //total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  //printf("Allocating %.5zu Gb for %d particles\n",total_memory,cp.npart);
  P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  assert(P != NULL);

  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  ind = 0;
  lee(filename,P,&ind);

  //fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);
}

extern void set_directory(const char * prefix, char * name, const type_int i, const type_int NNN, const type_real * fof)
{

  sprintf(name,"%s",snap.box);

  sprintf(name,"%s%.2d_%.1fMpc_%.2f_%.6d",name,Seg[i].flag,Seg[i].len/1000.0,Seg[i].razon,Seg[i].id);

  #ifdef NEW
    sprintf(name,"%s_new_%s",name,prefix);
  #else
    sprintf(name,"%s_%s",name,prefix);
  #endif

  #ifndef ORIGINAL
    sprintf(name,"%s_smooth",name);
  #endif

  sprintf(name,"%s_%.2d_%.4d",name,snap.num,NNN);

  #ifdef MCRITIC
	  sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%.2f_%.2f",name,fof[0],fof[1]);  

	return;
}

static void set_name(const char * prefix, char * name, const type_int NNN, const type_real * fof)
{
  sprintf(name,snap.seg);

  sprintf(name,"%s%.2d_%.4d",name,snap.num,NNN);

  sprintf(name,"%s_smooth",name);

  #ifdef NEW
    sprintf(name,"%s_new_%s",name,prefix);
  #else
    sprintf(name,"%s_%s",name,prefix);
  #endif

  #ifdef MCRITIC
	  sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);  

	return;
}

extern void read_segment(type_int NNN, type_real *fof)
{
  char  filename[200];
  FILE  *pf;
  type_int   i, j, k;
  type_real  r[3];
  type_real  cte_rvir;

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  GREEN("********** IMPORTANTE ***********\n");
  cp.simulbox  = cp.lbox*1e3;
  cp.simunpart = cp.npart;
  cte_rvir  = cbrt(3.0/(4.0*M_PI*cp.Mpart*(type_real)cp.simunpart));
  cte_rvir *= cp.simulbox;
  cte_rvir *= fof[1];
  fprintf(stdout,"lbox     %g Kpc\n",cp.simulbox);
  GREEN("**********************************\n");  

  set_name("segmentos",filename,NNN,fof);
  fprintf(stdout,"%s\n",filename);
  fflush(stdout);

  pf = fopen(filename,"rb"); 

  fread(&cp.nseg,sizeof(type_int),1,pf);

  fprintf(stdout,"Segmentos %d\n",cp.nseg);
  fflush(stdout);

  Seg = (struct segmentstd *) malloc(cp.nseg*sizeof(struct segmentstd));

  for(i=0;i<cp.nseg;i++)
  {
    fread(&Seg[i].size,sizeof(type_int),1,pf);
    
    Seg[i].Pos = (type_real *) malloc(3*Seg[i].size*sizeof(type_real));

    for(k=0;k<Seg[i].size;k++)
    {
      fread(r, sizeof(type_real), 3, pf);

      for(j=0;j<3;j++)
        Seg[i].Pos[3*k+j] = r[j];
    }
        
  }

  fclose(pf);

  set_name("propiedades",filename,NNN,fof);
  fprintf(stdout,"%s\n",filename);

  pf = fopen(filename,"rb"); 

  fread(&k,sizeof(type_int),1,pf);

  assert(k==cp.nseg);

  fprintf(stdout,"Propiedades Segmentos %d\n",cp.nseg);
  fflush(stdout);

  for(i=0;i<cp.nseg;i++)
  {  
    Seg[i].id = i;
    fread(&Seg[i].flag,sizeof(type_int),1,pf);
    fread(&k,sizeof(type_int),1,pf);
    fread(Seg[i].Mass,sizeof(type_real),2,pf);
    fread(Seg[i].Vnodos,sizeof(type_real),6,pf);
    fread(&Seg[i].razon,sizeof(type_real),1,pf);
    fread(&Seg[i].len,sizeof(type_real),1,pf);
    fread(&Seg[i].elong,sizeof(type_real),1,pf);
    fread(&Seg[i].rms,sizeof(type_real),1,pf);

    Seg[i].Rvir[0] = cbrt(Seg[i].Mass[0])*cte_rvir;
    Seg[i].Rvir[1] = cbrt(Seg[i].Mass[1])*cte_rvir;

    assert(k==Seg[i].size);
  }

  fclose(pf);

  j = 0;
  for(i=0;i<cp.nseg;i++)
  {
    if(Seg[i].flag == 0) continue;
    #ifdef CUT_IN_LEN
    if(Seg[i].len<LEN_MIN || Seg[i].len>LEN_MAX) continue;
    #endif

    Seg[j] = Seg[i];
    j++;
  }

  Seg = (struct segmentstd *) realloc(Seg,j*sizeof(struct segmentstd));
  cp.nseg = j;

  GREEN("********** IMPORTANTE ***********\n");
  sprintf(filename,"cut LEN %f Mpc REALOCATEA %d fil\n",0.5*(LEN_MIN+LEN_MAX)/1000.0f,cp.nseg);GREEN(filename);
  GREEN("**********************************\n");
  fflush(stdout);

  for(i=0;i<cp.nseg;i++)
  {

    if(i%200==0) 
    {
      fprintf(stdout,"%i %.4g\n",i,(type_real)i/(type_real)cp.nseg);
      fflush(stdout);
    }

    set_directory("box",filename,i,NNN,fof);

    sprintf(filename,"%s/box_gadget_%.6d.bin",filename,Seg[i].id);

    read_gadget(filename);

    // MUY IMPORTANTE
    change_positions();
  
    grid.nobj = cp.npart;
    grid.ngrid = (long)(cp.lbox/HSPH);

    if(grid.ngrid > NGRIDMAX)
    {
      grid.ngrid = NGRIDMAX;
    }
    
    grid_init();
    grid_build();

    init_sph(i);

    free(P);
    grid_free();

  }

  return;
}


