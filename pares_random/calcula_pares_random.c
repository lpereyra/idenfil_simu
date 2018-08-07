#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "colores.h"
#include "leesnap.h"
#include "calcula_pares_random.h"

#define NSIZE_MAX 1000


static int compare_descend(const void *a, const void *b)
{
  if(((struct sort_prop *) a)->mass > (((struct sort_prop *) b)->mass))
    return -1;

  if(((struct sort_prop *) a)->mass < (((struct sort_prop *) b)->mass))
    return +1;

  return 0;
}

static void reorder_grups(int *Id)
{
  int i;
  struct grup_data Grsave, Grsource;
  int idsource, idsave, dest;

  for(i=0;i<cp.ngrup;i++)
  {
      if(Id[i] != i)
	    {

    	  Grsource = Gr[i];
	      idsource = Id[i];
    	  dest = Id[i];

	      while(1)
	      {
	         Grsave = Gr[dest];
	         idsave = Id[dest];

	         Gr[dest] = Grsource;
	         Id[dest] = idsource;

	         if(dest == i)  break;

	         Grsource = Grsave;
	         idsource = idsave;

	         dest = idsource;

        }

	    }
   }

   return;
}

static void find_cut(int *k, const int NCut)
{
  const int Nitera = 10000;
  int a, b, i, j, m;

  a = 0;
  b = cp.ngrup;
  
  RED("Inicio Cut...\n");

  j = 0;

  while(j<Nitera)
  {
    m = (b+a)/2;

    if(Gr[m].NumPart>NCut)
    { 
      a = m;
    }else if(Gr[m].NumPart<NCut){
      b = m;
    }else{

      break; // sale del while
    }

    j++;
  }

  for(i=m; i<cp.ngrup; i++)
  {
    if(Gr[i].NumPart<NCut)
    {
      *k = i;
      break;
    }
  }

  GREEN("Finalizo Cut...\n");
  return;

}

static void sort_halos()
{
  int i;
  int *Id;
  struct sort_prop *mp;

  mp = (struct sort_prop *) malloc(sizeof(struct sort_prop) * cp.ngrup);
  Id = (int *) malloc(sizeof(int) * cp.ngrup);

  RED("Inicio Sort Masa...\n");
  fflush(stdout);

  for(i=0;i<cp.ngrup;i++)
  {
    mp[i].index = i;
    mp[i].mass  = cp.Mpart*Gr[i].NumPart;
  }

  qsort(mp, cp.ngrup, sizeof(struct sort_prop), compare_descend);

  for(i=0;i<cp.ngrup;i++)
    Id[mp[i].index] = i;

  free(mp);

  reorder_grups(Id);

  free(Id);

  GREEN("Finalizo Sort Masa...\n");

  return;
}

void crea_random(int NpartCut)
{
  int i,j,k,Npares,Tid,c,realloc_size;
  type_real Posprima[3];
  type_real lbox2 = cp.lbox*0.5;
  type_real a_mass, b_mass, dis;
  FILE *pfpares;
  char filename[200];
  struct data * data_share;

  k = c = 0;

  /// ordeno los halos por masa
  if(NpartCut!=0)
  {
    
    sort_halos();
    find_cut(&k,NpartCut);
  
  }else{

    // TODOS LOS GRUPOS
    k = cp.ngrup;
     
  }

  assert(k!=0);
  
  data_share = (struct data *) malloc(NTHREADS*sizeof(struct data));
  
  for(i=0;i<NTHREADS;i++)
  {
    data_share[i].nsize = 0;
    data_share[i].matrix = (struct info *) malloc(NSIZE_MAX*sizeof(struct info));
  }

  Npares = k*(k-1)/2;  

  sprintf(filename,"%.2d_%.4d_pares_halos_%.2f_%.2f.bin",snap.num,NpartCut,fof[0],fof[1]);
  pfpares = fopen(filename,"w");
  fwrite(&Npares,sizeof(int),1,pfpares);
  fprintf(stdout,"Total de Pares %d\n",Npares); 

  for(i=0;i<k-1;i++)
  { 

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(j,Tid, \
    Posprima,dis,realloc_size) shared(cp,i,k,Gr,data_share,lbox2) \
    reduction(-:Npares)
    for(j=i+1;j<k;j++)
    {

      Tid = omp_get_thread_num(); 

      Npares--;      
      Posprima[0] = Gr[i].Pos[0] - Gr[j].Pos[0];
      Posprima[1] = Gr[i].Pos[1] - Gr[j].Pos[1];
      Posprima[2] = Gr[i].Pos[2] - Gr[j].Pos[2];

      #ifdef PERIODIC
      if(Posprima[0] >  lbox2) Posprima[0] = Posprima[0] - cp.lbox;
      if(Posprima[1] >  lbox2) Posprima[1] = Posprima[1] - cp.lbox;
      if(Posprima[2] >  lbox2) Posprima[2] = Posprima[2] - cp.lbox;
      if(Posprima[0] < -lbox2) Posprima[0] = Posprima[0] + cp.lbox;
      if(Posprima[1] < -lbox2) Posprima[1] = Posprima[1] + cp.lbox;
      if(Posprima[2] < -lbox2) Posprima[2] = Posprima[2] + cp.lbox;
      #endif      

      dis    = sqrt(Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2]);

      if(dis<LEN_MIN || dis>LEN_MAX) continue;

      data_share[Tid].matrix[data_share[Tid].nsize].id_a = i;
      data_share[Tid].matrix[data_share[Tid].nsize].id_b = j;
      data_share[Tid].matrix[data_share[Tid].nsize].dist = dis;
      data_share[Tid].nsize++;

      if((data_share[Tid].nsize%NSIZE_MAX) == 0)
      {
        realloc_size = data_share[Tid].nsize/NSIZE_MAX;
        realloc_size *= NSIZE_MAX;
        realloc_size += NSIZE_MAX;

        data_share[Tid].matrix = (struct info *) \
                                 realloc(data_share[Tid].matrix, \
                                 realloc_size*sizeof(struct info));
      } 
    }
  }

  assert(Npares==0);

  for(Tid=0;Tid<NTHREADS;Tid++)
  { 
    c += data_share[Tid].nsize;

    for(k=0;k<data_share[Tid].nsize;k++)
    {
      i = data_share[Tid].matrix[k].id_a;
      j = data_share[Tid].matrix[k].id_b;

      a_mass = cp.Mpart*Gr[i].NumPart;
      b_mass = cp.Mpart*Gr[j].NumPart;
      dis    = data_share[Tid].matrix[k].dist;


      fwrite(&Gr[j].id,sizeof(int),1,pfpares);
      fwrite(&Gr[i].id,sizeof(int),1,pfpares);
      fwrite(&dis,sizeof(float),1,pfpares);
      fwrite(&b_mass,sizeof(float),1,pfpares);
      fwrite(&a_mass,sizeof(float),1,pfpares);
    }

    free(data_share[Tid].matrix);
  }

  rewind(pfpares);
  fwrite(&c,sizeof(int),1,pfpares);
  fclose(pfpares);

  fprintf(stdout,"Pares %d in MIN_MAX %f %f Mpc\n",c,LEN_MIN/1000.,LEN_MAX/1000.);

  return;

}

