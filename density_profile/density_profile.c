#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <mpi.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "colores.h"
#include "calcula.h"

#ifdef MARIO
  static void read_seg_mario(void);
#endif

#define DIV_CEIL(x,y) (x+y-1)/y

static void cut_reparte(int NumProc, int MyProc);


int main(int argc, char **argv)
{
  int NumProc, MyProc;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NumProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);
  double start,end;
#ifndef MARIO
  type_int NNN;
#endif
  TIMER(start);

  init_variables(argc,argv);
  omp_set_nested(1);
#ifndef MARIO
  NNN  = atoi(argv[2]);
#endif
 
  //Lee archivos de la simulacion //
  read_gadget();

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  fprintf(stdout,"lbox %g Kpc\n",cp.lbox);
  GREEN("**********************************\n");
#ifdef MARIO
  read_seg_mario();
#else
  read_segment(NNN,fof);
#endif

  limpia_calculados();
  cut_reparte(NumProc,MyProc);
  propiedades(fof);

  free(Seg);
  free(Gr);
  free(P);

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);

  MPI_Finalize();

  return(EXIT_SUCCESS);

}

#ifdef MARIO
static void read_seg_mario(void)
{
  type_int    i,j,k;
  type_real r[3];
  FILE *pf;
  char filename[200];
  
  sprintf(filename,"%s%s",snap.root_fil,snap.name_fil);
  fprintf(stdout,"Filename %s\n",filename);

  pf = fopen(filename,"rb"); 
  fread(&j,sizeof(type_int),1,pf);
  fprintf(stdout,"Segmentos %d\n",j);

  cp.ngrup=0;
  for(i=0;i<j;i++)
  {
    fread(&k,sizeof(type_int),1,pf);
    cp.ngrup+=k;
    fseek(pf,3*k*sizeof(type_real),SEEK_CUR);
  }  
  fclose(pf);

  Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));

  pf = fopen(filename,"rb"); 
  fread(&cp.nseg,sizeof(type_int),1,pf);
  fprintf(stdout,"Segmentos %d\n",cp.nseg);
  fflush(stdout);

  assert(j==cp.nseg);

  Seg = (struct segmentstd *) malloc(cp.nseg*sizeof(struct segmentstd));

  j=0;
  for(i=0;i<cp.nseg;i++)
  {
    fread(&Seg[i].size,sizeof(type_int),1,pf);
    
    Seg[i].id = i;
    Seg[i].start = j;
    Seg[i].Rvir_2[0] = -10.0;
    Seg[i].Rvir_2[1] = -10.0;

    for(k=0;k<Seg[i].size;k++)
    {
      fread(&r[0], sizeof(type_real), 3, pf);
      #ifdef SINTETICOS_RECTOS
        Gr[j].Pos[0] = r[0]*1000.0;
        Gr[j].Pos[1] = r[1]*1000.0;
        Gr[j].Pos[2] = r[2]*1000.0;
      #else
        Gr[j].Pos[0] = r[0];
        Gr[j].Pos[1] = r[1];
        Gr[j].Pos[2] = r[2];
      #endif
      j++;
    }
  }

  fclose(pf);

#ifdef CUT_ELONGACION

  float len, elongacion;
  type_int l, c;
  l = c = 0;

  for(i=0;i<cp.nseg;i++)
  {
    len = 0.0f;
    for(k=Seg[i].start+1;k<Seg[i].start+Seg[i].size;k++)
    {
      for(j=0;j<3;j++)
      {
        r[j] = Gr[k].Pos[j]-Gr[k-1].Pos[j];    
        #ifdef PERIODIC
        if(r[j]> 0.5*cp.lbox) r[j] -= cp.lbox;
        if(r[j]<-0.5*cp.lbox) r[j] += cp.lbox;
        #endif
      }
      len += sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    }

    for(j=0;j<3;j++)
    {
      r[j] = Gr[Seg[i].start+Seg[i].size-1].Pos[j]-Gr[Seg[i].start].Pos[j];    
      #ifdef PERIODIC
      if(r[j]> 0.5*cp.lbox) r[j] -= cp.lbox;
      if(r[j]<-0.5*cp.lbox) r[j] += cp.lbox;
      #endif
    }
    
    elongacion = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/len;
    
    if(elongacion<CUT_ELONG)
      continue;

    //if(Seg[i].size>80)
    //  fprintf(stdout,"-------- %d %d %f\n",i,Seg[i].size,len);

    for(k=Seg[i].start;k<Seg[i].start+Seg[i].size;k++)
    {
      Gr[c].Pos[0] = Gr[k].Pos[0];
      Gr[c].Pos[1] = Gr[k].Pos[1];
      Gr[c].Pos[2] = Gr[k].Pos[2];
      
      c++;
    }

    Seg[i].start = c-Seg[i].size;
    Seg[l] = Seg[i];

    l++;

  }

  cp.nseg = l;
  cp.ngrup = c;

  Seg = (struct segmentstd *) realloc(Seg,cp.nseg*sizeof(struct segmentstd));
  Gr  = (struct grup_data  *) realloc(Gr,cp.ngrup*sizeof(struct grup_data));

  fprintf(stdout,"Filamentos CUT ELONGATION %d\n",cp.nseg);
  fprintf(stdout,"Nodos      CUT ELONGATION %d\n",cp.ngrup);
  fprintf(stdout,"Segmentos  CUT ELONGATION %d\n",cp.ngrup-cp.nseg);

#endif

  return;
}

#endif

static void cut_reparte(int NumProc, int MyProc)
{
  type_int i,k,c,l;
  type_int start;
  type_int end;

  start = MyProc == 0         ? 0       :  MyProc   *DIV_CEIL(cp.nseg,NumProc);
  end   = MyProc == NumProc-1 ? cp.nseg : (MyProc+1)*DIV_CEIL(cp.nseg,NumProc);

  //c = MyProc == 0         ? 0                :  MyProc   *DIV_CEIL((cp.ngrup-cp.nseg),NumProc);
  //l = MyProc == NumProc-1 ? cp.ngrup-cp.nseg : (MyProc+1)*DIV_CEIL((cp.ngrup-cp.nseg),NumProc);
 
  //start = c;
  //end   = l;

  //for(i=0;i<cp.nseg;i++)
  //{
  //  if(c >= Seg[i].start-i)
  //    start = i;

  //  if(l > Seg[i].start-i)
  //    end = i;
  //}

  //start = MyProc == 0         ? 0       : start;
  //end   = MyProc == NumProc-1 ? cp.nseg : end;

  c = l = 0;
  for(i=start;i<end;i++)
  {
    for(k=Seg[i].start;k<Seg[i].start+Seg[i].size;k++)
    {
      Gr[c].Pos[0] = Gr[k].Pos[0];
      Gr[c].Pos[1] = Gr[k].Pos[1];
      Gr[c].Pos[2] = Gr[k].Pos[2];
      
      c++;
    }

    Seg[i].start = c-Seg[i].size;
    Seg[l] = Seg[i];

    l++;
  }

  cp.nseg = l;
  cp.ngrup = c;

  Seg = (struct segmentstd *) realloc(Seg,cp.nseg*sizeof(struct segmentstd));
  Gr  = (struct grup_data  *) realloc(Gr,cp.ngrup*sizeof(struct grup_data));

  fprintf(stdout,"Proc            %d\n",MyProc);
  fprintf(stdout,"Start           %d\n",start);
  fprintf(stdout,"End             %d\n",end);
  fprintf(stdout,"Filamentos CUT  %d\n",cp.nseg);
  fprintf(stdout,"Nodos      CUT  %d\n",cp.ngrup);
  fprintf(stdout,"Segmentos  CUT  %d\n",cp.ngrup-cp.nseg);
  fflush(stdout);

  return;

}
