#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

extern void init_variables(int argc, char **argv){
  FILE *pfin;
  char filename[200];
  type_int  i;

  RED("Initializing variables...\n");

  size_real = sizeof(type_real);
  fprintf(stdout,"type_real: %zu\n",size_real);

  size_int  = sizeof(type_int);
  fprintf(stdout,"type_int:  %zu\n",size_int);

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

  if(!fscanf(pfin,"%lf \n",&cp.soft))
  {
    sprintf(message,"can't read file `%s`\nneed softening of simulation\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%d  \n",&nfrac))
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

  if(!fscanf(pfin,"%d  \n",&nbins))
  {
    sprintf(message,"can't read file `%s`\nneed Nbins\n",filename);RED(message);
    exit(0);
  }

  #ifdef FIXED_SEPARATION

    #ifdef PRECDOUBLE
    if(!fscanf(pfin,"%lf  \n",&RLEN))
    #else
    if(!fscanf(pfin,"%f   \n",&RLEN))
    #endif
    {
      sprintf(message,"can't read file `%s`\n LEN MAX\n",filename);RED(message);
      exit(0);
    }
 
    #ifdef PRECDOUBLE
    if(!fscanf(pfin,"%lf  \n",&RSEP))
    #else
    if(!fscanf(pfin,"%f   \n",&RSEP))
    #endif
    {
      sprintf(message,"can't read file `%s`\n SEP MAX\n",filename);RED(message);
      exit(0);
    }

  #else

    if(!fscanf(pfin,"%d  \n",&ncil))
    {
      sprintf(message,"can't read file `%s`\nneed Ncil\n",filename);RED(message);
      exit(0);
    }

  #endif

  #ifdef MCRITIC

  #ifdef PRECDOUBLE
    if(!fscanf(pfin,"%lf \n",&m_critica))
  #else
    if(!fscanf(pfin,"%f \n",&m_critica))
  #endif
    {
      sprintf(message,"can't read file `%s`\nneed M_CRITIC\n",filename);RED(message);
      exit(0);
    }
  #endif

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
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"Snapnum:                 %d\n",snap.num);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Softening of simulation: %lf \n",cp.soft);BLUE(message);
  sprintf(message,"Identification steps:    %d\n",nfrac);BLUE(message);
  for(i=0;i<nfrac;i++)
  {
     sprintf(message,"%d overdensity          %.2f\n",i,fof[i]);BLUE(message);
     fof[i] = cbrt(1./(1.+fof[i]));
  }
  sprintf(message,"Num bins:               %d\n",nbins);BLUE(message);
  #ifdef FIXED_SEPARATION
    RED("FIXED SEPARATION\n");
    sprintf(message,"RLEN                %f Mpc\n",RLEN);BLUE(message);
    sprintf(message,"RSEP                %f Mpc\n",RSEP);BLUE(message);
  #else
    sprintf(message,"Num cylindres:          %d\n",ncil);BLUE(message);
  #endif

  #ifdef VEL_RELATIVA
    RED("CREA ARCHIVO PARA VELOCIDADES RELATIVAS\n");
  #endif
  #ifdef BIN_LOG
    RED("BIN LOG\n");
  #else
    RED("BIN LIN\n");
  #endif
  #ifdef EXTEND
    RED("EXTEND\n");
  #endif
  #ifdef REORDER
    RED("REORDER\n");
  #endif
  #ifdef COLUMN
    RED("COLUMN\n");
  #endif

  #ifdef MCRITIC
  sprintf(message,"M_CRIT [10^10 Msol / h]  %f\n",m_critica);RED(message);
  #endif

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

  GREEN("END\n");
}

#ifdef COLUMN

extern int allocate_particles(const type_int size)
{

  P.x  = (type_real *) malloc(size*sizeof(type_real));
  P.y  = (type_real *) malloc(size*sizeof(type_real));
  P.z  = (type_real *) malloc(size*sizeof(type_real));
  #ifdef STORE_VELOCITIES
  P.vx = (type_real *) malloc(size*sizeof(type_real));
  P.vy = (type_real *) malloc(size*sizeof(type_real));
  P.vz = (type_real *) malloc(size*sizeof(type_real));
  #endif
  #ifdef STORE_IDS
  P.id = (type_int  *) malloc(size*sizeof(type_int));
  #endif
  //P.gr = (type_int  *) malloc(size*sizeof(type_int));
  
  if(!P.x || !P.y || !P.z) 
  {
    fprintf(stderr, "cannot allocate pos particles\n" );
    return(0);
  }

  #ifdef STORE_VELOCITIES
  if(!P.vx || !P.vy || !P.vz) 
  {
    fprintf(stderr, "cannot allocate vel particles\n" );
    return(0);
  }
  #endif

  #ifdef STORE_IDS
  if(!P.id) 
  {
    fprintf(stderr, "cannot allocate ids particles\n" );
    return(0);
  }
  #endif

  return ( 1 );
}

extern void free_particles( void )
{
  if(P.x)  free(P.x);
  if(P.y)  free(P.y);
  if(P.z)  free(P.z);
  #ifdef STORE_VELOCITIES
  if(P.vx) free(P.vx);
  if(P.vy) free(P.vy);
  if(P.vz) free(P.vz);
  #endif   
  #ifdef STORE_IDS
  if(P.id) free(P.id);
  #endif   
}

#endif
