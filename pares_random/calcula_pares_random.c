#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "colores.h"
#include "leesnap.h"
#include "calcula_pares_random.h"

void find_cut(int *k, int NCut)
{
  int a, b, i, m;

  a = 0;
  b = cp.ngrup;
  
  while(1)
  {
    m = (b-a)/2 + a;

    if(Gr[m].NumPart>NCut)
    { 
      a = m;
    }else if(Gr[m].NumPart<NCut){
      b = m;
    }else{

      for(i=m; i<cp.ngrup; i++)
      {
        if(Gr[i].NumPart<NCut)
        {
          *k = i;
          break;
        }
      }

      break; // sale del while
    }

  }

  return;

}

void crea_random(int NpartCut)
{
  int i,j,k,c,Npares;
  type_real Posprima[3];
  type_real lbox2 = cp.lbox*0.5;
  type_real a_mass, b_mass, dis;
  FILE *pfpares;
  char filename[200];

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

  Npares = k*(k-1)/2;  

  sprintf(filename,"%.2d_%.4d_pares_halos_%.2f_%.2f.bin",snap.num,NpartCut,fof[0],fof[1]);
  pfpares = fopen(filename,"w");
  fwrite(&Npares,sizeof(int),1,pfpares);
  fprintf(stdout,"Total de Pares %d\n",Npares); 

  for(i=0;i<k-1;i++)
  { 

    a_mass = cp.Mpart*Gr[i].NumPart;

    for(j=i+1;j<k;j++)
    {
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

      if(dis<3000 || dis>15000) continue;

      b_mass = cp.Mpart*Gr[j].NumPart;

      fwrite(&Gr[j].id,sizeof(int),1,pfpares);
      fwrite(&Gr[i].id,sizeof(int),1,pfpares);
      fwrite(&dis,sizeof(float),1,pfpares);
      fwrite(&b_mass,sizeof(float),1,pfpares);
      fwrite(&a_mass,sizeof(float),1,pfpares);
      c++;
    }
  }

  rewind(pfpares);
  fwrite(&c,sizeof(int),1,pfpares);
  fclose(pfpares);

  fprintf(stdout,"Pares %d\n",c); 
  assert(Npares==0);

  return;

}

void sort_halos()
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

int compare_descend(const void *a, const void *b)
{
  if(((struct sort_prop *) a)->mass > (((struct sort_prop *) b)->mass))
    return -1;

  if(((struct sort_prop *) a)->mass < (((struct sort_prop *) b)->mass))
    return +1;

  return 0;
}

void reorder_grups(int *Id)
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

	      do
	      {
	         Grsave = Gr[dest];
	         idsave = Id[dest];

	         Gr[dest] = Grsource;
	         Id[dest] = idsource;

	         if(dest == i)  break;

	         Grsource = Grsave;
	         idsource = idsave;

	         dest = idsource;

        }while(1);

	    }
   }

   return;
}
