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
#include "prune.hh"

#ifdef BRANCH_SURVIVE
  type_int N_part_survive;
#endif

static void Write_Segments(const type_int NumPartCut, std::vector<std::pair<type_int,type_int> > &vec_orden, \
type_int * __restrict__ Padre, type_int * __restrict__ Rank, type_real * __restrict__ fof);

int main(int argc, char **argv)
{
  double start,end;
  double MassCut;
  type_int NumPartCut;
  type_int *Padre, *Rank;
  std::vector<std::vector<type_int> > adjacency_list;
  std::vector<std::pair<type_int,type_int> > vec_orden;

  TIMER(start);

  init_variables(argc,argv);
  omp_set_nested(1);

  // Lee archivos de la simulacion
  leeheader();

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

  read_grup_fof(fof[1]);

  read_mst(adjacency_list,fof);

	#ifdef PRUNED

	  Podado(LEVEL_PRUNED,adjacency_list);  

  	#ifdef BRANCH_SURVIVE
		  fprintf(stdout,"Sobreviven en ramas los nodos con %d\n",N_part_survive);
  		fflush(stdout);
  	#endif

  #endif
	
	Padre = (type_int *) malloc(cp.ngrup*sizeof(type_int));
  Rank =  (type_int *) malloc(cp.ngrup*sizeof(type_int));
  DL(vec_orden,adjacency_list,Padre,Rank);

  fprintf(stdout,"Escribe\n");
  fflush(stdout);

  Write_Segments(NumPartCut, vec_orden,Padre,Rank,fof);

  free(Gr);
  free(Padre);
  free(Rank);
  adjacency_list.clear();
  vec_orden.clear();

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}

static void set_name(const char * prefix, char * name, const type_int NumPartCut, const type_real * __restrict__ fof)
{

  sprintf(name,"%.2d_%.4d_%s",snap.num,NumPartCut,prefix);

  #ifdef MCRITIC
	  sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

	//#ifdef PRUNED
	//  sprintf(name,"%s_pruned_%.2d",name,LEVEL_PRUNED);
  //#endif

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);

	return;
}

static void Write_Segments(const type_int NumPartCut, std::vector<std::pair<type_int,type_int> > &vec_orden, \
type_int * __restrict__ Padre, type_int * __restrict__ Rank, type_real * __restrict__ fof)
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

	/*
  j = 0;
  while(!vec_orden.empty())
  {

    i = vec_orden.back().second;

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

			if(aux.size()!=1) // es igual a uno, no entra al while
			{
				segmentos.push_back(aux);
				j++;
			}

      aux.clear();
    }

    vec_orden.pop_back();
  }
	*/

  j = 0;
  while(!vec_orden.empty())
  {
    i = vec_orden.back().second;

    if(Padre[i]>=cp.ngrup || Rank[i] >= cp.ngrup)
		{
    	vec_orden.pop_back();
			continue;
		}

    aux.push_back(i);
    Rank[i] += cp.ngrup;
    id = Padre[i];

    while(id<cp.ngrup)
    {
      if(Rank[id]>=cp.ngrup)
      {
        aux.push_back(id);
    		Rank[id] += cp.ngrup;
        break;
      }

      if(Gr[id].NumPart>NumPartCut)
      {
        aux.push_back(id);
        Rank[id] += cp.ngrup;
        segmentos.push_back(aux);
        aux.clear();
        j++;

				if(Padre[id]>=cp.ngrup)
					break;
      }

      aux.push_back(id);
      Rank[id] += cp.ngrup;
      id = Padre[id];               
    }

		if(aux.size() > 1) // es igual a uno, no entra al while
		{
			segmentos.push_back(aux);
			j++;
		}

    aux.clear();
    vec_orden.pop_back();
  }

  assert(segmentos.size()==j);

  for(i=0;i<j;i++)
  {
	  if(Gr[segmentos[i][0]].NumPart>Gr[segmentos[i].back()].NumPart)
		{
	  	id = segmentos[i].size()-1;
			aux.push_back(segmentos[i][id]);
	  	for(k=segmentos[i].size()-1;k>1;k--)
	  	{
	  	  if(Gr[segmentos[i][id]].NumPart<Gr[segmentos[i][k-1]].NumPart)
				{
					aux.push_back(segmentos[i][k-1]);
					segmentos.push_back(aux);
					aux.clear();
					id = k-1;
					aux.push_back(segmentos[i][id]);
				}else{
					aux.push_back(segmentos[i][k-1]);
				}
	 		}
			aux.push_back(segmentos[i][0]);

			if(aux.size() != segmentos[i].size())
				segmentos[i].swap(aux);
			aux.clear();
			
		}else{

	  	id = 0;
			aux.push_back(segmentos[i][id]);
	  	for(k=1;k<segmentos[i].size()-1;k++)
	  	{
	  	  if(Gr[segmentos[i][id]].NumPart<Gr[segmentos[i][k]].NumPart)
				{
					aux.push_back(segmentos[i][k]);
					segmentos.push_back(aux);
					aux.clear();
					id = k;
					aux.push_back(segmentos[i][id]);
				}else{
					aux.push_back(segmentos[i][k]);
				}
	 		}
			aux.push_back(segmentos[i][segmentos[i].size()-1]);

			if(aux.size() != segmentos[i].size())
				segmentos[i].swap(aux);
			aux.clear();

		}
	}

	fprintf(stdout,"Seg Antes %u\n",j);
	fprintf(stdout,"Seg %lu\n",segmentos.size());
	j = segmentos.size();

	#ifdef NEW
		set_name("new_segmentos",filename,NumPartCut,fof);
	#else
		set_name("segmentos",filename,NumPartCut,fof);
	#endif
  pfout=fopen(filename,"w");
  fwrite(&j,sizeof(type_int),1,pfout);

	#ifdef NEW
		set_name("new_propiedades",filename,NumPartCut,fof);
	#else
		set_name("propiedades",filename,NumPartCut,fof);
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

    fwrite(&k,sizeof(type_int),1,pfpropiedades);

    k = (type_int)aux.size();
    fwrite(&k,sizeof(type_int),1,pfout);
    fwrite(&k,sizeof(type_int),1,pfpropiedades);

    fwrite(&Gr[aux[0]].Mass,sizeof(type_real),1,pfpropiedades); // ESCRIBE LA MASA DE LA PRIMERA PUNTA
    fwrite(&Gr[id].Mass,sizeof(type_real),1,pfpropiedades);     // ESCRIBE LA MASA DE LA SEGUNDA PUNTA

    fwrite(&Gr[aux[0]].vcm,sizeof(type_real),3,pfpropiedades);     // ESCRIBE LA VELOCIDAD DE LA SEGUNDA
    fwrite(&Gr[id].vcm,sizeof(type_real),3,pfpropiedades);     // ESCRIBE LA VELOCIDAD DE LA SEGUNDA

    r = Gr[aux[0]].Mass/Gr[id].Mass;
    fwrite(&r,sizeof(type_real),1,pfpropiedades); // ESCRIBE LA RAZON DE MASAS

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

    elong = sqrt(dux*dux+duy*duy+duz*duz);

    dux *= (1.0f/elong);
		duy *= (1.0f/elong);
		duz *= (1.0f/elong);

    lenr = rms = 0.0f;

    fwrite(&aux[0],sizeof(type_int),1,pfout);

    for(k=1;k<aux.size();k++)
    {
      fwrite(&aux[k],sizeof(type_int),1,pfout);

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
          
      r = dx*dux+dy*dux+dz*dux;    //dot
      r = dx*dx+dy*dy+dz*dz - r*r; //dist - dot
      r = fabs(r);                 //por errores de redondeo

      rms += r;
    }

    r = elong/lenr;

    k = (type_int)aux.size();
    rms /= (type_real)k;
    rms = sqrt(rms);

    fwrite(&lenr,sizeof(type_real),1,pfpropiedades);
    fwrite(&r,sizeof(type_real),1,pfpropiedades);
    fwrite(&rms,sizeof(type_real),1,pfpropiedades);

    aux.clear();
  }

  fclose(pfout);
  fclose(pfpropiedades);

  while(!segmentos.empty())
    segmentos.pop_back();

  fprintf(stdout,"Segmentos %d\n",j);

  return;

}
