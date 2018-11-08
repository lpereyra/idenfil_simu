#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
#include <vector>
#include "voro++.hh"

#include "colores.h"
#include "cosmoparam.h"
#include "variables.h"
#include "voronoi.h"
#include "colores.h"

static struct info *Pinfo; // Informacion de las Particulas
static type_real Delta[3]; // Grid en cada dimension
static type_real pmin[3];  // Min coordenadas
static type_real pmax[3];  // Max coordenadas

static void Distribution(const type_int Start, const type_int End, type_int * __restrict__ IsBad, type_int * __restrict__ IsBis)
{
  type_int i;
  type_int IndX,IndY,IndZ,SendProc;

  for(i=Start;i<End;i++)
  {

    IndX = (type_int) floorf(P[i].Pos[0]/Delta[0]); 
    IndY = (type_int) floorf(P[i].Pos[1]/Delta[1]); 
    IndZ = (type_int) floorf(P[i].Pos[2]/Delta[2]); 

  	SendProc = IndX + NX*(IndY + NY*IndZ);

    if(SendProc==NTHREADS) SendProc = NTHREADS-1;

    P[i].Save = -(SendProc+1);

    IsBis[i] = IsBad[i] = 1;

  }

  return;
}

static void Distribution_Shift(const type_int Tid, const type_int Start, const type_int End, type_int * __restrict__ IsBad, type_int * __restrict__ IsBis, const type_real Shift[])
{
  type_int i;
  type_int IndX,IndY,IndZ,SendProc;

  for(i=Start;i<End;i++)
  {
    if(P[i].Pos[0] >= cp.lbox) P[i].Pos[0] -= cp.lbox;     
    if(P[i].Pos[1] >= cp.lbox) P[i].Pos[1] -= cp.lbox;     
    if(P[i].Pos[2] >= cp.lbox) P[i].Pos[2] -= cp.lbox;    
    if(P[i].Pos[0] <  0.0)     P[i].Pos[0] += cp.lbox;     
    if(P[i].Pos[1] <  0.0)     P[i].Pos[1] += cp.lbox;     
    if(P[i].Pos[2] <  0.0)     P[i].Pos[2] += cp.lbox;    

    if(P[i].Save<0 && IsBad[i] == 0) P[i].Save *= -1;

    if(IsBis[i] == 0) continue;
    
    if(P[i].Pos[0] < fabs(Shift[0])*Delta[0]) P[i].Pos[0] += cp.lbox;
    if(P[i].Pos[1] < fabs(Shift[1])*Delta[1]) P[i].Pos[1] += cp.lbox;
    if(P[i].Pos[2] < fabs(Shift[2])*Delta[2]) P[i].Pos[2] += cp.lbox;

    IndX = (type_int) floorf(P[i].Pos[0]/Delta[0] - Shift[0]); 
    IndY = (type_int) floorf(P[i].Pos[1]/Delta[1] - Shift[1]); 
    IndZ = (type_int) floorf(P[i].Pos[2]/Delta[2] - Shift[2]); 

  	SendProc = IndX + NX*(IndY + NY*IndZ);

    if(SendProc==NTHREADS) SendProc = NTHREADS-1;

    if(P[i].Save>0)
      P[i].Save = SendProc+1;
    else
      P[i].Save = -(SendProc+1);

    #ifdef SAVENEIGH
      if(P[i].Save<0 && IsBad[i] == 1)
        P[i].neigh.clear();
    #endif
  }

  return;
}

static void Voronoi(const type_int Tid, const bool periodic, type_int * __restrict__ IsBad, type_int * __restrict__ IsBis, \
type_int * __restrict__ Malas, type_int * __restrict__ Duplicadas)
{
  
  const double gap = 1.0e-6;
  double       x_min,x_max,y_min,y_max,z_min,z_max;
  type_int     Npart,Ngrid,Flag;
  type_int     i,j,k,id,idv,idvv;
  std::vector<type_int>  vec,vvec;
  voro::voronoicell_neighbor cell;
  #ifdef SAVECENTROID
  std::vector<double> vertex;
  #endif

  if(!periodic)
  {

    // Busco limites de la caja
    x_min = y_min = z_min = 2.0*cp.lbox+gap;
    x_max = y_max = z_max = 0.0-gap;
     
    Npart = 0;
    for(i=0; i<cp.npart; i++) 
    {
      if(IsBis[i] == 0) continue;

      if(abs(P[i].Save)==Tid+1)
      {
        if(P[i].Pos[0] > x_max) x_max = P[i].Pos[0];
        if(P[i].Pos[1] > y_max) y_max = P[i].Pos[1];
        if(P[i].Pos[2] > z_max) z_max = P[i].Pos[2];
        if(P[i].Pos[0] < x_min) x_min = P[i].Pos[0];
        if(P[i].Pos[1] < y_min) y_min = P[i].Pos[1];
        if(P[i].Pos[2] < z_min) z_min = P[i].Pos[2];
        Npart++;
      }

    }

    x_min -= gap*(x_max-x_min);
    x_max += gap*(x_max-x_min);
    y_min -= gap*(y_max-y_min);
    y_max += gap*(y_max-y_min);
    z_min -= gap*(z_max-z_min);
    z_max += gap*(z_max-z_min);
   
  }else{

    Npart = 0;
    for(i=0; i<cp.npart; i++) 
    {
      if(IsBis[i] == 0) continue;

      if(abs(P[i].Save)==Tid+1)
        Npart++;
    }

    x_min = y_min = z_min = 0.0;	  
    x_max = y_max = z_max = cp.lbox;	  

  }

  if(Npart==0)
  {
    RED("Warning Npart = 0!!!\n");
    return;
  }
  
  // Creo el voro::container
  Ngrid = (type_int) pow((double)Npart/5.0,1.0/3.0);
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,\
                      Ngrid,Ngrid,Ngrid,periodic,periodic,periodic,8);

  for(i=0; i<cp.npart; i++) 
  {
    if(IsBis[i] == 0) continue;

    if(abs(P[i].Save)==Tid+1)
    {
      con.put(i,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);	  
      IsBis[i] = IsBad[i] = 0;
    }
  }

  assert(Npart==con.total_particles());

  fprintf(stdout,"Proc %d: %d NumPart (xbool = %s, ybool = %s, zbool = %s)\n",Tid,Npart,\
  periodic ? "true" : "false", periodic ? "true" : "false", periodic ? "true" : "false");
  fflush(stdout);

  // Creo clase loop y guardo informacion del grid y orden en el grid
  voro::c_loop_all clo(con);

  if(clo.start()) do {
     id = clo.pid();
     Pinfo[id].ijk = clo.ijk;
     Pinfo[id].q = clo.q;
  } while (clo.inc());
  // Guardo informacion voronoi y localizo particulas a recalcular

  *Malas = *Duplicadas = 0;
  if(clo.start()) do {
   
    id = clo.pid();

    if(P[id].Save>0) continue;

    if(!con.compute_cell(cell,clo)) continue;

    cell.neighbors(vec);   
    Flag = 0;

    for(k=0; k<(type_int)vec.size(); k++) 
    {
       if(vec[k] < 0)
       {
         Flag++;
         break;
       }
    }
  
    if(Flag != 0)
    {

      if(IsBad[id] == 0)
      {
        IsBad[id] = 1;
        (*Malas)++;
      }

      if(IsBis[id] == 0)
      {
        IsBis[id] = 1;
        (*Duplicadas)++;
      }
      
      // Vecinos  
      for(k=0; k<(type_int)vec.size(); k++) 
      {

        if(vec[k] < 0) continue;

        idv = vec[k];    

        if(IsBad[idv] == 0)
        {
          IsBad[idv] = 1;
          (*Malas)++;
        }

        if(IsBis[idv] == 0)
        {
          IsBis[idv] = 1;
          (*Duplicadas)++;
        }
         
        // Vecinos de vecinos
        if(!con.compute_cell(cell,Pinfo[idv].ijk,Pinfo[idv].q)) continue;
      
        cell.neighbors(vvec);

        for(j=0; j<(type_int)vvec.size(); j++) 
        {

          if(vvec[j] < 0)
          {
            Flag++;
            continue;
          }

          idvv = vvec[j];

          if(IsBis[idvv] == 0)
          {
            IsBis[idvv] = 1;
            (*Duplicadas)++;
          }

        }

        vvec.clear(); 

      }

     }else{

       if(IsBad[id] == 0)
       {	  

         P[id].vol = cell.volume();

         #ifdef SAVENEIGH
           P[id].neigh.swap(vec);
         #endif

         #ifdef SAVECENTROID
          cell.vertices(vertex);
       
          P[id].c[0] = P[id].c[1] = P[id].c[2] = 0.0f;

          for(i=0;i<cell.p;i++)
          {  
            P[id].c[0] += vertex[3*i]/(double)cell.p;
            P[id].c[1] += vertex[3*i+1]/(double)cell.p;
            P[id].c[2] += vertex[3*i+2]/(double)cell.p;
          }

         #endif
       }

     }

     vec.clear(); 

  }while (clo.inc());

  con.clear();

  return;

}

static void Voronoi_End(const type_int Tid, const bool periodic, type_int * __restrict__ IsBad, type_int * __restrict__ IsBis, \
type_int * __restrict__ Malas, type_int * __restrict__ Duplicadas)
{
  
  const double gap = 1.0e-6;
  double       x_min,x_max,y_min,y_max,z_min,z_max;
  type_int     Npart,Ngrid,Flag;
  type_int     i,j,k,id,idv,idvv;
  voro::voronoicell_neighbor cell;
  #ifdef SAVENEIGH
  std::vector<type_int>  vec;
  #endif
  #ifdef SAVECENTROID
  std::vector<double> vertex;
  #endif

  if(!periodic)
  {

    // Busco limites de la caja
    x_min = y_min = z_min = 2.0*cp.lbox+gap;
    x_max = y_max = z_max = 0.0-gap;
     
    Npart = 0;
    for(i=0; i<cp.npart; i++) 
    {
      if(IsBis[i] == 0) continue;

      if(abs(P[i].Save)==Tid+1)
      {
        if(P[i].Pos[0] > x_max) x_max = P[i].Pos[0];
        if(P[i].Pos[1] > y_max) y_max = P[i].Pos[1];
        if(P[i].Pos[2] > z_max) z_max = P[i].Pos[2];
        if(P[i].Pos[0] < x_min) x_min = P[i].Pos[0];
        if(P[i].Pos[1] < y_min) y_min = P[i].Pos[1];
        if(P[i].Pos[2] < z_min) z_min = P[i].Pos[2];
        Npart++;
      }

    }

    x_min -= gap*(x_max-x_min);
    x_max += gap*(x_max-x_min);
    y_min -= gap*(y_max-y_min);
    y_max += gap*(y_max-y_min);
    z_min -= gap*(z_max-z_min);
    z_max += gap*(z_max-z_min);
   
  }else{

    Npart = 0;
    for(i=0; i<cp.npart; i++) 
    {
      if(IsBis[i] == 0) continue;

      if(abs(P[i].Save)==Tid+1)
        Npart++;
    }

    x_min = y_min = z_min = 0.0;	  
    x_max = y_max = z_max = cp.lbox;	  

  }

  if(Npart==0)
  {
    RED("Warning Npart = 0!!!\n");
    return;
  }
  
  // Creo el voro::container
  Ngrid = (type_int) pow((double)Npart/5.0,1.0/3.0);
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,\
                      Ngrid,Ngrid,Ngrid,periodic,periodic,periodic,8);

  for(i=0; i<cp.npart; i++) 
  {
    if(IsBis[i] == 0) continue;

    if(abs(P[i].Save)==Tid+1)
    {
      con.put(i,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);	  
      IsBis[i] = IsBad[i] = 0;
    }
  }

  assert(Npart==con.total_particles());

  fprintf(stdout,"Proc %d: %d NumPart (xbool = %s, ybool = %s, zbool = %s)\n",Tid,Npart,\
  periodic ? "true" : "false", periodic ? "true" : "false", periodic ? "true" : "false");
  fflush(stdout);

  // Creo clase loop y guardo informacion del grid y orden en el grid
  voro::c_loop_all clo(con);

  if(clo.start()) do {
   
    id = clo.pid();

    if(P[id].Save>0) continue;

    if(!con.compute_cell(cell,clo)) continue;

    P[id].vol = cell.volume();

    #ifdef SAVENEIGH
      cell.neighbors(vec);   
      P[id].neigh.swap(vec);
     vec.clear(); 
    #endif

    #ifdef SAVECENTROID
     cell.vertices(vertex);
    
     P[id].c[0] = P[id].c[1] = P[id].c[2] = 0.0f;

     for(i=0;i<cell.p;i++)
     {  
       P[id].c[0] += vertex[3*i]/(double)cell.p;
       P[id].c[1] += vertex[3*i+1]/(double)cell.p;
       P[id].c[2] += vertex[3*i+2]/(double)cell.p;
     }

    #endif

  }while (clo.inc());

  con.clear();

  return;

}

static void change_positions()
{
  type_int i;

  for(i=0;i<3;i++)
  {
    pmin[i] = 1.E26; 
    pmax[i] = -1.E26;
  }

  for(i=0;i<cp.npart;i++)
  {
    if(P[i].Pos[0] > pmax[0]) pmax[0] = P[i].Pos[0];
    if(P[i].Pos[0] < pmin[0]) pmin[0] = P[i].Pos[0];
    if(P[i].Pos[1] > pmax[1]) pmax[1] = P[i].Pos[1];
    if(P[i].Pos[1] < pmin[1]) pmin[1] = P[i].Pos[1];
    if(P[i].Pos[2] > pmax[2]) pmax[2] = P[i].Pos[2];
    if(P[i].Pos[2] < pmin[2]) pmin[2] = P[i].Pos[2];
  }

  fprintf(stdout,"xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  fprintf(stdout,"ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  fprintf(stdout,"zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);
 
  for(i=0;i<cp.npart;i++)
  {
    P[i].Pos[0] -= pmin[0];
    P[i].Pos[1] -= pmin[1];
    P[i].Pos[2] -= pmin[2];
  }

  cp.lbox = 0.0f;
  for(i = 0; i < 3; i++)
    if(cp.lbox < (pmax[i] - pmin[i])) cp.lbox = (pmax[i] - pmin[i]);

  cp.lbox *= 1.001;

  return;

}

static void rechange_positions()
{
  type_int i;

  for(i=0;i<cp.npart;i++)
  {
    #ifdef LLOYD
      P[i].Pos[0] += P[i].c[0]+pmin[0];
      P[i].Pos[1] += P[i].c[1]+pmin[1];
      P[i].Pos[2] += P[i].c[2]+pmin[2];
    #else
      P[i].Pos[0] += pmin[0];
      P[i].Pos[1] += pmin[1];
      P[i].Pos[2] += pmin[2];
    #endif

    #ifdef PERIODIC
      if(P[i].Pos[0] >= cp.lbox) P[i].Pos[0] -= cp.lbox;     
      if(P[i].Pos[1] >= cp.lbox) P[i].Pos[1] -= cp.lbox;     
      if(P[i].Pos[2] >= cp.lbox) P[i].Pos[2] -= cp.lbox;    
      if(P[i].Pos[0] <  0.0)     P[i].Pos[0] += cp.lbox;     
      if(P[i].Pos[1] <  0.0)     P[i].Pos[1] += cp.lbox;     
      if(P[i].Pos[2] <  0.0)     P[i].Pos[2] += cp.lbox;    
    #endif
  }

  return;

}

extern void init_voro(void)
{
  type_int i,j,k,idim,Tid,Flag;
  type_int *Start, *End;
  type_int *Malas, *Duplicadas;
  type_int *IsBad, *IsBis;
  type_real Shift[3];         // El desplazamiento en cada dimension

  // ES MUY IMPORTANTE
  change_positions();

  Delta[0] = cp.lbox/(type_real)NX;	
  Delta[1] = cp.lbox/(type_real)NY;	
  Delta[2] = cp.lbox/(type_real)NZ;	

  Pinfo       = (struct info *) malloc(cp.npart*sizeof(struct info));
  IsBad       = (type_int *) malloc(cp.npart*sizeof(type_int));
  IsBis       = (type_int *) malloc(cp.npart*sizeof(type_int));
  Start       = (type_int *) malloc(NTHREADS*sizeof(type_int));
  End         = (type_int *) malloc(NTHREADS*sizeof(type_int));
  Malas       = (type_int *) malloc(NTHREADS*sizeof(type_int));
  Duplicadas  = (type_int *) malloc(NTHREADS*sizeof(type_int));
  
  Flag = cp.npart;

  RED("Abre la region paralela\n");
  sprintf(message,"%d THREADS\n",NTHREADS);
  sprintf(message,"computing %d cell\n",Flag);
  RED(message);
  fflush(stdout);

  #pragma omp parallel num_threads(NTHREADS) default(none) \
  private(i,j,k,idim,Tid) shared(stdout,cp,Flag,Start,End,IsBad,\
  IsBis,Malas,Duplicadas,Delta,Shift,P)
  {

    Tid = omp_get_thread_num();

    Start[Tid] = Tid==0 ? 0 : Tid*floor((type_real)cp.npart/NTHREADS);
    End[Tid] = Tid==NTHREADS-1 ? cp.npart : (Tid+1)*floor((type_real)cp.npart/NTHREADS);
    Malas[Tid] = Duplicadas[Tid] = 0;

    if(Flag != 0)
      Distribution(Start[Tid],End[Tid],IsBis,IsBad);

    #pragma omp barrier

    if(Tid==0)
    {
      j = 0;
      for(i=0;i<NTHREADS;i++)
        j+=End[i]-Start[i];

      if(j!=cp.npart)
      {
        fprintf(stdout,"\n%d NumPart != %d cp.npart\n",j,cp.npart);
        fflush(stdout);
        exit(0);
      }

      BLUE("Termina de Distribuir las Particulas\n");
      GREEN("Comienza Voronoi\n");
      fflush(stdout);

      Flag = j;
    }

    #pragma omp barrier

    if(Flag!=0)
      Voronoi(Tid,false,IsBad,IsBis,&Malas[Tid],&Duplicadas[Tid]);

    #pragma omp barrier

    if(Tid==0)
    {
      j = k = 0;
      for(i=0;i<NTHREADS;i++) 
      {
        k += Malas[i];
        j += Duplicadas[i];
      }

      Flag = j;

      fprintf(stdout,"Recalcular = %d ---> %f %% | Necesarias = %d ---> %f %% \n",
      k,100.0*((type_real)k/(type_real)cp.npart),j,100.0*((type_real)j/(type_real)cp.npart)); 
    }  
    
    #ifdef PERIODIC

    #pragma omp barrier

    for(idim=0;idim<4;idim++)
    {

      if(Tid==0)
      {
        Shift[0] = Shift[1] = Shift[2] = 0.0;

        if(idim==3)
          Shift[0] = Shift[1] = Shift[2] = 0.5;
        else
          Shift[idim] = 0.5;

        BLUE("Comienza a Distribuir las Particulas Shift\n");
        fprintf(stdout,"Shift = (%f,%f,%f)\n",Shift[0],Shift[1],Shift[2]);
        fflush(stdout);
      }

      #pragma omp barrier

      if(Flag != 0)
        Distribution_Shift(Tid,Start[Tid],End[Tid],IsBad,IsBis,Shift);

      #pragma omp barrier

      if(Tid==0)
      {
        BLUE("Termina de Distribuir las Particulas\n");
        GREEN("Comienza Voronoir Shift\n");
        fflush(stdout);
      }

      #pragma omp barrier
 
      if(Flag!=0)
        Voronoi(Tid,false,IsBad,IsBis,&Malas[Tid],&Duplicadas[Tid]);

      #pragma omp barrier

      if(Tid==0)
      {
        j = k = 0;
        for(i=0;i<NTHREADS;i++) 
        {
          k += Malas[i];
          j += Duplicadas[i];
        }

        fprintf(stdout,"Recalcular = %d ---> %f %% | Necesarias = %d ---> %f %% \n",
        k,100.0*((type_real)k/(type_real)cp.npart), 
        j,100.0*((type_real)j/(type_real)cp.npart)); 
        fflush(stdout);

        Flag = j;
      }
      
    }

    #endif

    for(i=Start[Tid];i<End[Tid];i++)
    {
      if(P[i].Pos[0] >= cp.lbox) P[i].Pos[0] -= cp.lbox;     
      if(P[i].Pos[1] >= cp.lbox) P[i].Pos[1] -= cp.lbox;     
      if(P[i].Pos[2] >= cp.lbox) P[i].Pos[2] -= cp.lbox;   
      if(P[i].Pos[0] <  0.0)     P[i].Pos[0] += cp.lbox;     
      if(P[i].Pos[1] <  0.0)     P[i].Pos[1] += cp.lbox;     
      if(P[i].Pos[2] <  0.0)     P[i].Pos[2] += cp.lbox;    
 
      if(P[i].Save<0 && IsBad[i] == 0) P[i].Save *= -1;

      if(IsBis[i] == 0) continue;

      if(P[i].Save > 0) 
        P[i].Save =  1;
      else
        P[i].Save = -1;

      #ifdef SAVENEIGH
        if(P[i].Save<0 && IsBad[i] == 1)
          P[i].neigh.clear();
      #endif
    }


  }

  RED("Cierra la region paralela\n");
  fflush(stdout);
 
  if(Flag!=0)
  {
    GREEN("Voronoi Serial sobre las particulas que quedan\n");
    fflush(stdout);
    #ifdef PERIODIC
      Voronoi_End(0,true,IsBad,IsBis,&Malas[0],&Duplicadas[0]);
    #else
      Voronoi_End(0,false,IsBad,IsBis,&Malas[0],&Duplicadas[0]);
    #endif
  }

  free(IsBad);
  free(IsBis);
  free(Start);
  free(End);
  free(Malas);
  free(Duplicadas);
  free(Pinfo);

  rechange_positions();

  GREEN("Finaliza Voronoi\n");
  fflush(stdout);

  return;

}
