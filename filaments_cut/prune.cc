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
#include "cosmoparam.hh"
#include "leesnap.hh"
#include "prune.hh"

#ifdef PRUNED

	static type_int DFS(const type_int i, std::vector<std::vector<type_int> > &vec, type_int cut)
	{
	  type_int j,k,id,idv;
	
	  j = 1;
	  id = i;
	  idv = vec[id][0];
	
	  while(vec[idv].size()==2)
	  {
	
	    #ifdef BRANCH_SURVIVE
	      if(Gr[idv].NumPart>N_part_survive)
	      {
	        j += cut+1;
	      }
	    #endif
	   
	    k   = vec[idv][0]==id ? vec[idv][1] : vec[idv][0];
	    id  = idv;
	    idv = k;
	    j++;        
	  }
		
		#ifdef BRANCH_SURVIVE
	    if(Gr[idv].NumPart>N_part_survive)
	    {
	      j += cut+1;
	    }
	  #endif
	
	  if(j<=cut)
	    return 1;
	  else 
	    return 0;
	
	}
	
	static void Delete_Branch(const type_int i, std::vector<std::vector<type_int> > &vec)
	{
	  type_int k, id, idv;
	
	  id = i;
	  idv = vec[id][0];
	  vec[id].clear();
	 
	  while(vec[idv].size()==2)
	  {
	
	    k = vec[idv][0]==id ? vec[idv][1] : vec[idv][0];
	
	    vec[idv].clear();
	
	    id = idv;
	    idv = k;
	  }

	  if(vec[idv].size()==1)
	    vec[idv].clear();

	  return;
	}
	
	extern void Podado(const type_int level, std::vector<std::vector<type_int> > &vec)
	{
	  type_int i, k, cut, N_threads, itera;
	  std::vector<type_int> aux;
	
	  #ifdef NTHREADS
	  omp_set_dynamic(0);
	  omp_set_num_threads(NTHREADS);
	  #endif
	
	  N_threads = NTHREADS;
	
	  fprintf(stdout,"level podado %d\n",level);
	  fflush(stdout);
	
	  cut = 1;
	  //cut = level;
	
	  do{
	
	    itera = 0;
	
	    #ifdef BRANCH_SURVIVE
	      #pragma omp parallel for num_threads(N_threads) schedule(static) default(none) \
	      shared(cp,Gr,vec,cut,itera,N_part_survive) 
	    #else
	      #pragma omp parallel for num_threads(N_threads) schedule(static) default(none) \
	      shared(cp,Gr,vec,cut,itera)  
	    #endif
	    for(i=0;i<cp.ngrup;i++)
	    {
	
	      if(vec[i].size()==1)
	      {
	        
	      #ifdef BRANCH_SURVIVE
	        if(Gr[i].NumPart>N_part_survive) continue;
	      #endif
	
	        if(DFS(i,vec,cut))       
	        {
	          Delete_Branch(i,vec);
	          itera = 1;
	        }
	  
	      }
	
	    }
	
	    if(itera==0) break;
	
	    #pragma omp parallel for num_threads(N_threads) schedule(static) default(none) \
	    private(aux,k) shared(cp,vec) 
	    for(i=0;i<cp.ngrup;i++)
	    {
	      if(vec[i].size()>2)
	      {
	        for(k=0;k<vec[i].size();k++)
	        {
	          if(!vec[vec[i][k]].empty())
	            aux.push_back(vec[i][k]);
	        }
	
	        swap(aux,vec[i]);
	        aux.clear();
	
	      }
	    }
	
	    if(cut==level) break;
	
	    cut++;
	
	  }while(1);
	
	  fprintf(stdout,"Termina el podado\n");
	  fflush(stdout);
	
	  return;
	  
	}

#endif

#ifdef NEW

	static type_real Calcula_Weight(const type_int v, const type_int u)
	{
		type_real w;
		type_real vdir[3];

	  vdir[0] = Gr[v].Pos[0] - Gr[u].Pos[0];
  	vdir[1] = Gr[v].Pos[1] - Gr[u].Pos[1];
    vdir[2] = Gr[v].Pos[2] - Gr[u].Pos[2];

	  #ifdef PERIODIC
  	vdir[0] = vdir[0] >= cp.lbox*0.5 ? vdir[0]-cp.lbox : vdir[0];
    vdir[0] = vdir[0] < -cp.lbox*0.5 ? vdir[0]+cp.lbox : vdir[0];

	  vdir[1] = vdir[1] >= cp.lbox*0.5 ? vdir[1]-cp.lbox : vdir[1];
  	vdir[1] = vdir[1] < -cp.lbox*0.5 ? vdir[1]+cp.lbox : vdir[1];

	  vdir[2] = vdir[2] >= cp.lbox*0.5 ? vdir[2]-cp.lbox : vdir[2];
  	vdir[2] = vdir[2] < -cp.lbox*0.5 ? vdir[2]+cp.lbox : vdir[2];
    #endif

    w = -((type_real)Gr[u].NumPart*(type_real)Gr[v].NumPart)/(vdir[0]*vdir[0]+vdir[1]*vdir[1]+vdir[2]*vdir[2]);

		return w; 

	}

#endif


#ifdef NEW
	static void DLU(const type_int id, const type_int pre, std::vector<std::vector<type_int> > &adj, type_int *Padre, type_int *Rank, type_real * Weight)
#else
	static void DLU(const type_int id, const type_int pre, std::vector<std::vector<type_int> > &adj, type_int *Padre, type_int *Rank)
#endif
{
  type_int k, idv;

  Padre[id] = pre;
  Rank[id]  = adj[id].size();

  for(k=0;k<Rank[id];k++)
  {
    idv = adj[id][k];
    if(idv != pre)
	  {
			#ifdef NEW
				Weight[idv] = Weight[id] + Calcula_Weight(id,idv);
      	DLU(idv,id,adj,Padre,Rank,Weight);
			#else
      	DLU(idv,id,adj,Padre,Rank);
			#endif
	  }
  }

  adj[id].clear();

  return;

}

extern void DL(std::vector<std::pair<type_int,type_int> > &vec_sort, std::vector<std::vector<type_int> > &vec, \
type_int * __restrict__ Padre, type_int * __restrict__ Rank)
{
  type_int i;
	#ifdef NEW
	type_real *Weight;
	#endif

  for(i=0;i<cp.ngrup;i++)
  {
    Padre[i] = cp.ngrup;
    Rank[i]  = 0;

    if(vec[i].size()>0)
      vec_sort.push_back(std::make_pair(Gr[i].NumPart,i));
  }

  sort(vec_sort.begin(),vec_sort.end());
	#ifdef NEW
	Weight = (type_real *) calloc(cp.ngrup,sizeof(type_real));
 	#endif

  for(i=vec_sort.size();i>0;i--)
  {
    type_int k = vec_sort[i-1].second;

    if(Padre[k]==cp.ngrup)
		{
			#ifdef NEW
      DLU(k,Padre[k],vec,Padre,Rank,Weight);
			#else
      DLU(k,Padre[k],vec,Padre,Rank);
			#endif
  	}
  }

	#ifdef NEW
  while(!vec_sort.empty())
    vec_sort.pop_back();

  for(i=0;i<cp.ngrup;i++)
  {
    if(Rank[i]>0)
      vec_sort.push_back(std::make_pair(Weight[i],i)); 
  }

  sort(vec_sort.begin(),vec_sort.end(), std::greater<std::pair<type_real, type_int> > ());

	free(Weight);
	#endif

  return;
}
