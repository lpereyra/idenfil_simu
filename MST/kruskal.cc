#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include <vector>
#include <functional>
#include <algorithm>

#include "variables.hh"
#include "leesnap.hh"
#include "kruskal.hh"

static type_int Root(type_int i, type_int *Padre)
{

 if(i != Padre[i])
   Padre[i] = Root(Padre[i],Padre);

 return Padre[i];
}

extern void Kruskal(type_int * __restrict__ Padre, type_int * __restrict__ Rank, \
std::vector<std::pair<type_real,std::pair<type_int,type_int> > > &edges)
{
  type_int j,k,id,idv,count;
  char filename[200];
  FILE *pfout;

  #ifdef MCRITIC
  sprintf(filename,"%.2d_mst_cut_%.2f_%.2f_%.2f.bin",snap.num,m_critica,fof[0],fof[1]);
  #else
  sprintf(filename,"%.2d_mst_%.2f_%.2f.bin",snap.num,fof[0],fof[1]);
  #endif

  count = 0;
  pfout=fopen(filename,"w");
  fwrite(&count,sizeof(type_int),1,pfout);

  sort(edges.begin(),edges.end(), std::greater<std::pair<type_real,std::pair<type_int,type_int> > >());

  while(!edges.empty()) 
  {
    j = edges.back().second.first;
    k = edges.back().second.second;
    edges.pop_back();

    id = Root(j,Padre);
    idv = Root(k,Padre);

    if(id!=idv)
    {
      if(Rank[id] < Rank[idv])
        Padre[id] = idv;
      else if(Rank[idv] < Rank[id])
        Padre[idv] = id;
      else
      {
        Padre[idv] = id;
        Rank[id]++;
      }

      fwrite(&j,sizeof(type_int),1,pfout);
      fwrite(&k,sizeof(type_int),1,pfout);
      count++;
    } 

  }

  fprintf(stdout,"%u NumEdges MST\n",count);
  fflush(stdout);

  rewind(pfout);
  fwrite(&count,sizeof(type_int),1,pfout);
  fclose(pfout);

  return;
}
