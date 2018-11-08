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

	fprintf(stdout,"Min %f  ",edges[0].first);
	fprintf(stdout,"Max %f\n",edges[edges.size()-1].first);
	fflush(stdout);

	type_real min_edges   = edges[0].first;
	type_real max_edges   = edges[edges.size()-1].first;
	type_real range_edges = max_edges-min_edges;

	for(j=0;j<edges.size();j++)
	{
		edges[j].first = (((double)edges[j].first-min_edges)/(double)range_edges);
	}

  //sort(edges.begin(),edges.end(), std::greater<std::pair<type_real,std::pair<type_int,type_int> > >());
  sort(edges.begin(),edges.end());

	fprintf(stdout,"Min %f  ",edges[0].first);
	fprintf(stdout,"Max %f\n",edges[edges.size()-1].first);
	fflush(stdout);

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
