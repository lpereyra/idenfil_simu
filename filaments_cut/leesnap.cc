#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <vector>
#include <assert.h>
#include "variables.hh"
#include "cosmoparam.hh"
#include "colores.hh"
#include "leesnap.hh"

static struct io_header header;
struct particle_data *P;
type_int *Index;
struct grup_data *Gr;

extern void leeheader(void)
{
  FILE *pf;
  char filename[200];
  type_int d1,d2;

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

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(filename,"Reescala\n");GREEN(filename);
  GREEN("**********************************\n");

	return;
}

extern void read_grup_fof(type_real prefix)
{
  char filename[200];
  type_int i;
  FILE *pfin;

  RED("Lectura de los grupos\n");

  #ifdef MCRITIC
    sprintf(filename,"%.2d_%.2f_centros_cut_%.2f.bin",snap.num,prefix,m_critica);
  #else
  	sprintf(snap.root,"../");
    sprintf(filename,"%s%.2d_%.2f_centros.bin",snap.root,snap.num,prefix);
  #endif

  pfin = fopen(filename,"rb"); 
  fread(&cp.ngrup,sizeof(type_int),1,pfin);
  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

  for(i=0;i<cp.ngrup;i++)
  {
    fread(&Gr[i].save,sizeof(type_int),1,pfin);
    fread(&Gr[i].id,sizeof(type_int),1,pfin);
    fread(&Gr[i].Pos[0],sizeof(type_real),1,pfin);
    fread(&Gr[i].Pos[1],sizeof(type_real),1,pfin);
    fread(&Gr[i].Pos[2],sizeof(type_real),1,pfin);
    fread(&Gr[i].NumPart,sizeof(type_int),1,pfin);
  }

  #ifdef MCRITIC
    fprintf(stdout,"Grupos %u mayores a mcrit %g\n",cp.ngrup,m_critica*1.e10);   
  #else
    fprintf(stdout,"Grupos %u\n",cp.ngrup); 
  #endif

  RED("Finaliza la lectura de los grupos\n"); fflush(stdout);

  fclose(pfin);

}


extern void read_mst(std::vector<std::vector<type_int> > &adjacency_list, type_real * __restrict__ fof)
{
  char filename[200];
  type_int i,j,k,len;
  FILE *pfin;

  RED("Allocating Adjacency List..\n");

  for(i=0;i<cp.ngrup;i++)
    adjacency_list.push_back(std::vector<type_int>());

  RED("Read MST\n");

  #ifdef MCRITIC
    sprintf(filename,"%.2d_mst_cut_%.2f_%.2f_%.2f.bin",snap.num,m_critica,fof[0],fof[1]);
  #else
    sprintf(filename,"%.2d_mst_%.2f_%.2f.bin",snap.num,fof[0],fof[1]);
  #endif

  pfin = fopen(filename,"rb"); 
  assert(pfin != NULL);

  fread(&len,sizeof(type_int),1,pfin);

  for(i=0;i<len;i++)
  {
    fread(&j,sizeof(type_int),1,pfin);
    fread(&k,sizeof(type_int),1,pfin);
    adjacency_list[j].push_back(k);
    adjacency_list[k].push_back(j);
  }

  RED("End Read MST\n"); fflush(stdout);

  fclose(pfin);

}
