#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "variables.hh"
#include "cosmoparam.hh"
#include "timer.hh"
#include "colores.hh"
#include "leesnap.hh"
#include "voronoi.hh"
#include "kruskal.hh"

int main(int argc, char **argv)
{
  type_int  i;
  type_int *Padre, *Rank;
  double start,end;
  std::vector<std::pair<type_real,std::pair<type_int,type_int> > > edges;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  //Lee archivos de la simulacion
  read_gadget();

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
  }

  Kruskal(Padre,Rank,edges);

  free(Gr);
  free(Padre);
  free(Rank);
  edges.clear();

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}
