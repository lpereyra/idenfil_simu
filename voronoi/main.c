#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "timer.h"
#include "colores.h"
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "voronoi.h"

static void Write();

int main(int argc, char **argv)
{
  double start,end;
  #ifdef LLOYD
  int i;

    #ifndef SAVECENTROID
      RED("Error ACTIVATE SAVECENTROID!!!!!!!!!\n");
      exit(0);
    #endif

  #endif

  TIMER(start);
  
  init_variables(argc,argv);

  if(NX*NY*NZ!=NTHREADS)
  {
    RED("Error!!!!!!!!!\n");
    RED("NX*NY*NZ != NTHREADS\n");
    fprintf(stdout,"%d != %d\n",NX*NY*NZ,NTHREADS);
    exit(0);
  }

  // Lee archivos de la simulacion 
  read_gadget();

  BLUE("Run...\n");
  fflush(stdout);

  #ifdef LLOYD
  for(i=0;i<100;i++)
  #endif
  init_voro();

  BLUE("Comienza a Escribir\n");
  fflush(stdout);
  Write();
  BLUE("Termina de Escribir\n");
  fflush(stdout);

  free(P);
 
  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}

static void Write()
{
  type_int i,id;
  const double volexa = cp.lbox*cp.lbox*cp.lbox;
  double voltot = 0.0f; 
  char filename[200];
  FILE *pfout;
  #ifdef WRITE_ASCII
   FILE *pfout_ascii;
  #endif
  #ifdef SAVENEIGH
    type_int j,k;
    FILE *pfout_neigh;
    #ifdef WRITE_ASCII
     FILE *pfout_neigh_ascii;
    #endif
  #endif
  #ifdef SAVECENTROID
    FILE *pfout_centroid;
    #ifdef WRITE_ASCII
     FILE *pfout_centroid_ascii;
    #endif
  #endif
 
  sprintf(filename,"volume.bin");
  pfout = fopen(filename,"w");
  fwrite(&cp.npart,sizeof(type_int),1,pfout);          
  #ifdef WRITE_ASCII
    sprintf(filename,"volume.dat");
    pfout_ascii = fopen(filename,"w");
    fprintf(pfout_ascii,"%d\n",cp.npart);         
  #endif

  #ifdef SAVECENTROID
  sprintf(filename,"centroid.bin");
  pfout_centroid = fopen(filename,"w");
  fwrite(&cp.npart,sizeof(type_int),1,pfout_centroid);
  #ifdef WRITE_ASCII
    sprintf(filename,"centroid.dat");
    pfout_centroid_ascii = fopen(filename,"w");
    fprintf(pfout_centroid_ascii,"%d\n",cp.npart);         
  #endif
  #endif

  #ifdef SAVENEIGH
  sprintf(filename,"vecinos.bin");
  pfout_neigh = fopen(filename,"w");
  fwrite(&cp.npart,sizeof(type_int),1,pfout_neigh);
  #ifdef WRITE_ASCII
    sprintf(filename,"vecinos.dat");
    pfout_neigh_ascii = fopen(filename,"w");
    fprintf(pfout_neigh_ascii,"%d\n",cp.npart);         
  #endif
  #endif

  for(i=0;i<cp.npart;i++)
  {
    #ifdef STORE_IDS
      id = P[i].id;
    #else
      id = i;
    #endif

    fwrite(&id,sizeof(type_int),1,pfout);          
    fwrite(&P[i].vol,sizeof(double),1,pfout);         

    #ifdef WRITE_ASCII
      fprintf(pfout_ascii,"%d %lf\n",id,P[i].vol);         
    #endif

    #ifdef SAVECENTROID

      P[i].c[0] += P[i].Pos[0];
      P[i].c[1] += P[i].Pos[1];
      P[i].c[2] += P[i].Pos[2];
 
      #ifdef PERIODIC
        if(P[i].c[0] >= cp.lbox) P[i].c[0] -= cp.lbox;     
        if(P[i].c[1] >= cp.lbox) P[i].c[1] -= cp.lbox;     
        if(P[i].c[2] >= cp.lbox) P[i].c[2] -= cp.lbox;    
        if(P[i].c[0] <  0.0)     P[i].c[0] += cp.lbox;     
        if(P[i].c[1] <  0.0)     P[i].c[1] += cp.lbox;     
        if(P[i].c[2] <  0.0)     P[i].c[2] += cp.lbox;    
      #endif

      fwrite(&id,sizeof(type_int),1,pfout_centroid);         
      fwrite(P[i].c,sizeof(double),3,pfout_centroid);        

      #ifdef WRITE_ASCII
        fprintf(pfout_centroid_ascii,"%d %lf %lf %lf\n",id,P[i].c[0],P[i].c[1],P[i].c[2]);         
      #endif

    #endif
      
    #ifdef SAVENEIGH

      k = (type_int)P[i].neigh.size();          
      fwrite(&id,sizeof(type_int),1,pfout_neigh);         
      fwrite(&k,sizeof(type_int),1,pfout_neigh);

      assert(k!=0);

      #ifdef WRITE_ASCII
        fprintf(pfout_neigh_ascii,"%d %d",id,k);         
      #endif

      for(j=0;j<k;j++)
      {
        #ifdef STORE_IDS
        fwrite(&P[P[i].neigh[j]].id,sizeof(type_int),1,pfout_neigh);          
        #ifdef WRITE_ASCII
          fprintf(pfout_neigh_ascii," %d",P[P[i].neigh[j]].id);         
        #endif
        #else
        fwrite(&P[i].neigh[j],sizeof(type_int),1,pfout_neigh);          
        #ifdef WRITE_ASCII
          fprintf(pfout_neigh_ascii," %d",P[i].neigh[j]);         
        #endif
        #endif
      }

      #ifdef WRITE_ASCII
       fprintf(pfout_neigh_ascii,"\n");         
      #endif
     
      P[i].neigh.clear();

    #endif

    voltot += P[i].vol;
  }

  fprintf(stdout,"vol real  %g\n",volexa);
  fprintf(stdout,"vol aprox %g\n",voltot);
  fprintf(stdout,"diff      %g\n",voltot-volexa);
  fflush(stdout);

  fclose(pfout);
  #ifdef WRITE_ASCII
   fclose(pfout_ascii);
  #endif

  #ifdef SAVENEIGH
    fclose(pfout_neigh);
    #ifdef WRITE_ASCII
      fclose(pfout_neigh_ascii);
    #endif
  #endif

  #ifdef SAVECENTROID
    fclose(pfout_centroid);
    #ifdef WRITE_ASCII
      fclose(pfout_centroid_ascii);
    #endif
  #endif
  
  return;
}
