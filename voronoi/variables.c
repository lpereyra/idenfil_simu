#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "voronoi.h"
#include "leesnap.h"
#include "cosmoparam.h"
#include "colores.h"

char   message[200]    ;
struct cosmoparam   cp ;
struct SnapST      snap;
struct particle_data *P;

extern void init_variables(int argc, char **argv)
{
  FILE *pfin;
  char filename[200];

  RED("Initializing variables...\n");

  sprintf(filename,"%s",argv[1]);
  if(!(pfin=fopen(filename,"r")))
  {
    sprintf(message,"can't open file `%s` \n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%d  \n",&snap.nfiles))
  {
    sprintf(message,"can't read file `%s`\nneed # of snapshots\n",filename);RED(message);
    exit(0);
  }
    
  if(!fscanf(pfin,"%s  \n",snap.root))
  {
    sprintf(message,"can't read file `%s`\nneed snapshots directory\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%s  \n",snap.name))
  {
    sprintf(message,"can't read file `%s`\nneed snapname\n",filename);RED(message);
    exit(0);
  }

  fclose(pfin);
  
  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  BLUE("********** Makefile Options ***********\n");
  #ifdef PERIODIC
  BLUE("  PERIODIC\n");
  #endif
  #ifdef PRECDOUBLE
  BLUE("  PRECDOUBLE\n");
  #endif
  #ifdef LONGIDS
  BLUE("  LONGIDS\n");
  #endif
  #ifdef POSFACTOR
  sprintf(message,"  POSFACTOR = %f\n",POSFACTOR);BLUE(message);
  #endif
  #ifdef VELFACTOR
  sprintf(message,"  VELFACTOR = %f\n",VELFACTOR);BLUE(message);
  #endif
  #ifdef STORE_VELOCITIES
  BLUE("  STORE_VELOCITIES\n");
  #endif
  #ifdef STORE_IDS
  BLUE("  STORE_IDS\n");
  #endif
  #ifdef SAVENEIGH
  BLUE("  SAVENEIGH\n");
  #endif
  #ifdef SAVECENTROID
  BLUE("  SAVECENTROID\n");
  #endif
  #ifdef WRITE_ASCII
  BLUE("  WRITE_ASCII\n");
  #endif
  sprintf(message,"  NX = %d\n",NX);BLUE(message);
  sprintf(message,"  NY = %d\n",NY);BLUE(message);
  sprintf(message,"  NZ = %d\n",NZ);BLUE(message);
  sprintf(message,"  NTHREADS = %d\n",NTHREADS);BLUE(message);
  GREEN("END\n");

}
