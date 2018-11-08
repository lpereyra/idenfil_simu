#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "variables.h"
#include "leesnap.h"
#include "cosmoparam.h"
#include "colores.h"
#include "grid.h"

char   message[200]    ;
struct cosmoparam   cp ;
struct SnapST      snap;
struct particle_data *P;
struct segmentstd  *Seg;
struct gridst      grid;
type_int          nfrac;
type_real          *fof;

extern void init_variables(int argc, char **argv)
{
  type_int i;
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

  if(!fscanf(pfin,"%u  \n",&nfrac))
  {
    sprintf(message,"can't read file `%s`\nneed identification steps\n",filename);RED(message);
    exit(0);
  }

  fof = (type_real *) malloc(nfrac*sizeof(type_real));

  for(i=0;i<nfrac;i++)
  {
    #ifdef PRECDOUBLE
    if(!fscanf(pfin,"%lf  \n",&fof[i]))
    #else
    if(!fscanf(pfin,"%f   \n",&fof[i]))
    #endif
    {
      sprintf(message,"can't read file `%s`\nneed %d identification step\n",filename,i);RED(message);
      exit(0);
    }
  }

  if(!fscanf(pfin,"%s  \n",snap.box))
  {
    sprintf(message,"can't read file `%s`\nneed box snapshots directory\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%s  \n",snap.seg))
  {
    sprintf(message,"can't read file `%s`\nneed seg snapname\n",filename);RED(message);
    exit(0);
  }

  fclose(pfin);

  /////////////////////////////////////////////////////////////////////////
  char *p = snap.name;
  while (*p) 
  { 
    if(isdigit(*p)) { // Upon finding a digit, ...
        snap.num = strtol(p, &p, 10); // Read a number, ...
    } else { // Otherwise, move on to the next character.
        p++;
    }
  }
  ///////////////////////////////////////////////////////////////////////

  BLUE("********** Information ***********\n");
  sprintf(message,"Snapnum:                 %d\n",snap.num);BLUE(message);
  sprintf(message,"Snapshots BOX directory: %s\n",snap.box);BLUE(message);
  sprintf(message,"Snapshots SEG directory: %s\n",snap.seg);BLUE(message);
  BLUE("********** Makefile Options ***********\n");
  for(i=0;i<nfrac;i++)
  {
     sprintf(message,"%d overdensity          %.2f\n",i,fof[i]);BLUE(message);
     fof[i] = cbrt(1./(1.+fof[i]));
  }
 
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
  #ifdef HSPH
  sprintf(message,"  H_SPH = %f Kpc\n",HSPH);BLUE(message);
  #endif
  #ifdef STORE_VELOCITIES
  BLUE("  STORE_VELOCITIES\n");
  #endif
  #ifdef STORE_IDS
  BLUE("  STORE_IDS\n");
  #endif
  #ifdef WRITE_ASCII
  BLUE("  WRITE_ASCII\n");
  #endif
  GREEN("END\n");
  fflush(stdout);

  return;

}
