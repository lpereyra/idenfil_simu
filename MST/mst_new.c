#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "voronoi.h"
#include "mst_kruskal.h"
#include "colores.h"

type_int  NumPartCut;
std::vector<std::pair<type_int,type_int> > mass_orden;
#ifdef BRANCH_SURVIVE
  type_int N_part_survive;
#endif

void Write_Segments(type_int *Padre, type_int *Rank, type_real *fof);

int main(int argc, char **argv)
{
  type_int  i;
  type_int  *Padre, *Rank;
  double start,end;
  double MassCut;
  std::vector<std::vector<type_int> > adjacency_list;
  std::vector<std::pair<type_real,std::pair<type_int,type_int> > > edges;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  /*Lee archivos de la simulacion*/
  read_gadget();

  MassCut = atof(argv[2]);
  NumPartCut = (type_int)(MassCut/cp.Mpart) ;  // Masa de la partícula [10^10 Msol / h]
  #ifdef BRANCH_SURVIVE
  N_part_survive = NumPartCut;
  #endif
  
  GREEN("********** Important *************\n");
  sprintf(message,"Mpart %g\n",cp.Mpart*1.e10);RED(message);
  sprintf(message,"Mass cut %g -> Nodos Npart %d Mass %g\n", \
  MassCut*1.e10,NumPartCut,NumPartCut*cp.Mpart*1e10);RED(message);
  GREEN("**********************************\n");
  fflush(stdout);

  select_particles_fof(fof[0]);

  read_grup_fof(fof[1]);

  Voronoi_Grupos(fof[0],edges);

  fprintf(stdout,"%lu NumEdges\n",edges.size());
  fflush(stdout);

  Padre = (type_int *) malloc(cp.ngrup*sizeof(type_int));
  Rank =  (type_int *) malloc(cp.ngrup*sizeof(type_int));

  for(i=0;i<cp.ngrup;i++)
  {
    Padre[i] = i;
    Rank[i] = 0;
    adjacency_list.push_back(std::vector<type_int>());
  }

  Kruskal(Padre,Rank,edges,adjacency_list);

  Podado(4,adjacency_list);  

  #ifdef BRANCH_SURVIVE
  fprintf(stdout,"Sobreviven en ramas los nodos con %d\n",N_part_survive);
  fflush(stdout);
  #endif

  for(i=0;i<cp.ngrup;i++)
  {
    Padre[i] = cp.ngrup;
    Rank[i]  = 0;

    if(adjacency_list[i].size()>0)
      mass_orden.push_back(std::make_pair(Gr[i].NumPart,i));
  }

  sort(mass_orden.begin(),mass_orden.end());
  
  for(i=mass_orden.size();i>0;i--)
  {
    type_int k = mass_orden[i-1].second;

    if(Padre[k]==cp.ngrup)
      DLU(k,Padre[k],adjacency_list,Padre,Rank);
  }

  fprintf(stdout,"Escribe\n");
  fflush(stdout);

  Write_Segments(Padre,Rank,fof);

  free(Gr);
  free(Padre);
  free(Rank);
  adjacency_list.clear();

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}

void Write_Segments(type_int *Padre, type_int *Rank, type_real *fof)
{
  char filename[200];
  FILE *pfout, *pfpropiedades;
  type_int i,j,k,id;
  type_real dx,dy,dz;
  type_real dux,duy,duz;
  type_real rms,elong;
  type_real r,lenr;
  std::vector<type_int> aux;
  std::vector<std::vector<type_int> > segmentos;

  j = 0;
  while(!mass_orden.empty())
  {

    i = mass_orden.back().second;

    //if(Padre[i]>=0)
    //if(Rank[i]==1 && Padre[i]>=0)
    if(Padre[i]<cp.ngrup)
    {
      aux.push_back(i);
      id = Padre[i];
      Rank[i] += cp.ngrup;

      while(id<cp.ngrup)
      {
        if(Rank[id]>=cp.ngrup)
        {
          aux.push_back(id);
          break;
        }

        //if(Rank[id]>2)
        if(Gr[id].NumPart>NumPartCut && Padre[id]!=cp.ngrup)
        {
          aux.push_back(id);

          segmentos.push_back(aux);
  
          aux.clear();
          j++;
        }

        aux.push_back(id);
        Rank[id] += cp.ngrup;
        id = Padre[id];               
      }

      segmentos.push_back(aux);

      aux.clear();
      j++;
    }

    mass_orden.pop_back();
  }

  assert(segmentos.size()==j);

  #ifdef MCRITIC
  sprintf(filename,"%.2d_%.4d_segmentos_cut_%.2f_%.2f_%.2f.bin",snap.num,NumPartCut,m_critica,fof[0],fof[1]);
  #else
  sprintf(filename,"%.2d_%.4d_segmentos_%.2f_%.2f.bin",snap.num,NumPartCut,fof[0],fof[1]);
  #endif
  pfout=fopen(filename,"w");
  fwrite(&j,sizeof(type_int),1,pfout);

  #ifdef MCRITIC
  sprintf(filename,"%.2d_%.4d_propiedades_cut_%.2f_%.2f_%.2f.bin",snap.num,NumPartCut,m_critica,fof[0],fof[1]);
  #else
  sprintf(filename,"%.2d_%.4d_propiedades_%.2f_%.2f.bin",snap.num,NumPartCut,fof[0],fof[1]);
  #endif

  pfpropiedades=fopen(filename,"w");
  fwrite(&j,sizeof(type_int),1,pfpropiedades);

  #ifdef SORT_DERECHA

  BLUE("********** Importante ***********\n");
  sprintf(filename,"Ordena el nodo de la derecha es más grande\n");RED(filename);
  BLUE("*********************************\n");

  #endif

  for(i=0;i<j;i++)
  {
    #ifdef SORT_DERECHA
    if(Gr[segmentos[i][0]].NumPart>Gr[segmentos[i].back()].NumPart)
      reverse(segmentos[i].begin(),segmentos[i].end());
    #endif

    aux = segmentos[i];

    id = aux.back();

    k = 0;
    if(Gr[aux[0]].NumPart>NumPartCut) k++;
    if(Gr[id].NumPart>NumPartCut)     k++;

    fwrite(&k,sizeof(int),1,pfpropiedades);

    k = (int)aux.size();
    fwrite(&k,sizeof(int),1,pfout);
    fwrite(&k,sizeof(int),1,pfpropiedades);

    r = (float)Gr[aux[0]].NumPart/(float)Gr[id].NumPart;
    fwrite(&r,sizeof(float),1,pfpropiedades); // ESCRIBE LA RAZON DE MASAS

    dux = Gr[id].Pos[0] - Gr[aux[0]].Pos[0];
    duy = Gr[id].Pos[1] - Gr[aux[0]].Pos[1];
    duz = Gr[id].Pos[2] - Gr[aux[0]].Pos[2];

    #ifdef PERIODIC
    dux = dux >= cp.lbox*0.5 ? dux-cp.lbox : dux;
    dux = dux < -cp.lbox*0.5 ? dux+cp.lbox : dux;

    duy = duy >= cp.lbox*0.5 ? duy-cp.lbox : duy;
    duy = duy < -cp.lbox*0.5 ? duy+cp.lbox : duy;

    duz = duz >= cp.lbox*0.5 ? duz-cp.lbox : duz;
    duz = duz < -cp.lbox*0.5 ? duz+cp.lbox : duz;
    #endif

    elong = dux*dux+duy*duy+duz*duz;

    lenr = rms = 0.0;

    fwrite(&aux[0],sizeof(int),1,pfout);

    for(k=1;k<aux.size();k++)
    {
      fwrite(&aux[k],sizeof(int),1,pfout);

      dx = Gr[aux[k]].Pos[0] - Gr[aux[k-1]].Pos[0];
      dy = Gr[aux[k]].Pos[1] - Gr[aux[k-1]].Pos[1];
      dz = Gr[aux[k]].Pos[2] - Gr[aux[k-1]].Pos[2];

      #ifdef PERIODIC
      dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
      dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

      dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
      dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

      dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
      dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
      #endif
      
      r = sqrt(dx*dx+dy*dy+dz*dz);

      lenr += r;

      if(k==aux.size()-1) continue;

      dx = Gr[aux[k]].Pos[0] - Gr[aux[0]].Pos[0];
      dy = Gr[aux[k]].Pos[1] - Gr[aux[0]].Pos[1];
      dz = Gr[aux[k]].Pos[2] - Gr[aux[0]].Pos[2];

      #ifdef PERIODIC
      dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
      dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

      dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
      dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

      dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
      dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
      #endif
          
      r = pow(dy*duz-dz*duy,2);
      r += pow(dz*dux-dx*duz,2);
      r += pow(dx*duy-dy*dux,2);
      r /= elong;
      
      rms += r;

    }

    r = sqrt(elong)/lenr;

    assert(r<1.+1e-06);

    k = (int)aux.size();
    rms /= (float)k;
    rms = sqrt(rms);

    fwrite(&lenr,sizeof(float),1,pfpropiedades);
    fwrite(&r,sizeof(float),1,pfpropiedades);
    fwrite(&rms,sizeof(float),1,pfpropiedades);

    aux.clear();
  }

  fclose(pfout);
  fclose(pfpropiedades);

  while(!segmentos.empty())
    segmentos.pop_back();

  fprintf(stdout,"Segmentos %d\n",j);

  return;

}
