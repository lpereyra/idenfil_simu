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

void read_segment_random(int NNN, type_real *fof)
{
  char  filename[200];
  int   flag,save,size,N0,N1;
  int   i,j,k,l,n;
  type_real leng;
  FILE  *pf;

  cp.nseg = 0;
 
  for(i=0;i<NTHREADS;i++)
  {
    #ifdef MCRITIC
    sprintf(filename,"../random_segment/%.2d_%.4d_random_fil_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NNN,m_critica,fof[0],fof[1],i);
    #else
    sprintf(filename,"../random_segment/%.2d_%.4d_random_fil_%.2f_%.2f.%.2d.bin",snap.num,NNN,fof[0],fof[1],i);
    #endif

    pf = fopen(filename,"rb"); 

    fread(&k,sizeof(int),1,pf);
    fread(&n,sizeof(int),1,pf);

    if(i==0) 
      cp.nrand = n;
    else
      assert(cp.nrand==n);

    cp.nseg += k;
   
    fprintf(stdout,"Tid %d Nfil %d\n",i,k);  

    fclose(pf);
  }

  PropSeg = (struct propsegmentstd *) malloc(cp.nseg*sizeof(struct propsegmentstd));
  Segrand = (struct segmentrandstd *) malloc(cp.nrand*cp.nseg*sizeof(struct segmentrandstd));

  for(i=0;i<NTHREADS;i++)
  {
    #ifdef MCRITIC
    sprintf(filename,"../random_segment/%.2d_%.4d_random_fil_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NNN,m_critica,fof[0],fof[1],i);
    #else
    sprintf(filename,"../random_segment/%.2d_%.4d_random_fil_%.2f_%.2f.%.2d.bin",snap.num,NNN,fof[0],fof[1],i);
    #endif

    pf = fopen(filename,"rb"); 

    fread(&k,sizeof(int),1,pf);
    fread(&n,sizeof(int),1,pf);

    for(j=0;j<k;j++)
    {
      fread(&save,sizeof(int),1,pf);
      fread(&flag,sizeof(int),1,pf);
      fread(&size,sizeof(int),1,pf);
      fread(&leng,sizeof(type_real),1,pf);
      fread(&N0,sizeof(int),1,pf);
      fread(&N1,sizeof(int),1,pf);

      PropSeg[save].flag = flag;
      PropSeg[save].size = size;
      PropSeg[save].len  = leng;
      PropSeg[save].M0   = (type_real)N0*cp.Mpart;
      PropSeg[save].M1   = (type_real)N1*cp.Mpart;

      for(l=0;l<n;l++)
      {
        Segrand[n*save+l].save = save;
        Segrand[n*save+l].Pos  = (type_real *) malloc(3*PropSeg[save].size*sizeof(type_real));

        fread(&Segrand[n*save+l].rand_ang[0],sizeof(type_real),3,pf);
        fread(&Segrand[n*save+l].Pos[0],sizeof(type_real),3*PropSeg[save].size,pf);
      }

    }
    
    fclose(pf);
  }

  return;
}
