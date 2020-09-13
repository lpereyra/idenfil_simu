#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "propiedades.h"
#include "colores.h"

extern void init_variables(int argc, char **argv)
{
  FILE *pfin;
  char filename[200];
  int i;

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

  if(!fscanf(pfin,"%lf \n",&cp.soft))
  {
    sprintf(message,"can't read file `%s`\nneed softening of simulation\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%u  \n",&nfrac))
  {
    sprintf(message,"can't read file `%s`\nneed identification steps\n",filename);RED(message);
    exit(0);
  }

  fof = (type_real *) malloc(nfrac*sizeof(type_real));

#ifdef READ_LEVELS

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

#else

  for(i = 0; i < nfrac; i++)
  {
    fof[i]  = 1.0f/(type_real)(nfrac + 1);
    fof[i] *= (type_real)(nfrac + 1 - i);
    fof[i] *= cbrt(1./(1.+FOF_OVERDENSITY));
  }

#endif

  fclose(pfin);

  /////////////////////////////////////////////////////////////////////////
  char *p = snap.name;
  while (*p) 
  { 
    if(isdigit(*p)) 
    { // Upon finding a digit, ...
        snap.num = strtol(p, &p, 10); // Read a number, ...
    } else { // Otherwise, move on to the next character.
        p++;
    }
  }
  ///////////////////////////////////////////////////////////////////////

  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"Snapname Num:            %d\n",snap.num);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Softening of simulation: %lf \n",cp.soft);BLUE(message);
  sprintf(message,"Identification steps:    %d\n",nfrac);BLUE(message);
#ifdef READ_LEVELS
  BLUE("************* Read Levels ************\n");
  for(i=0;i<nfrac;i++)
  {
    sprintf(message,"%d overdensity %.2f\n",i,fof[i]);BLUE(message);
    fof[i] = cbrt(1./(1.+fof[i]));
  }
#endif

  BLUE("********** Makefile Options ***********\n");
  #ifdef DEBUG
  BLUE("  DEBUG\n");
  #endif
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
  #ifdef LOCK
  BLUE("  USING LOCKS\n");
  #endif

  GREEN("END\n");
}

extern void write_properties(FILE *pfout, struct propiedades_st Prop)
{
	type_int   jj;

  fwrite(&Prop.npart,sizeof(type_int),1,pfout);
  fwrite(&Prop.pcm,sizeof(type_real),3,pfout);
#ifdef STORE_VELOCITIES  
  fwrite(&Prop.vcm,sizeof(type_real),3,pfout);
  fwrite(&Prop.mostbound,sizeof(type_real),3,pfout);
  fwrite(&Prop.sig,sizeof(type_real),3,pfout);
  fwrite(&Prop.L,sizeof(type_real),3,pfout);
  fwrite(&Prop.lambda,sizeof(type_real),1,pfout);
  fwrite(&Prop.m200,sizeof(type_real),1,pfout);
  fwrite(&Prop.r200,sizeof(type_real),1,pfout);
  fwrite(&Prop.v200,sizeof(type_real),1,pfout);
  fwrite(&Prop.mvir,sizeof(type_real),1,pfout);
  fwrite(&Prop.rvir,sizeof(type_real),1,pfout);
  fwrite(&Prop.vvir,sizeof(type_real),1,pfout);
  fwrite(&Prop.vmax,sizeof(type_real),1,pfout);
  fwrite(&Prop.Ep,sizeof(type_real),1,pfout);
  fwrite(&Prop.Ec,sizeof(type_real),1,pfout);
#endif  
  fwrite(&Prop.aa,sizeof(type_real),1,pfout);
  fwrite(&Prop.bb,sizeof(type_real),1,pfout);
  fwrite(&Prop.cc,sizeof(type_real),1,pfout);
  for(jj=0;jj<3;jj++)
    fwrite(&Prop.evec[jj],sizeof(type_real),3,pfout);
#ifdef STORE_VELOCITIES  
  fwrite(&Prop.aa_vel,sizeof(type_real),1,pfout);
  fwrite(&Prop.bb_vel,sizeof(type_real),1,pfout);
  fwrite(&Prop.cc_vel,sizeof(type_real),1,pfout);
  for(jj=0;jj<3;jj++)
    fwrite(&Prop.evec_vel[jj],sizeof(type_real),3,pfout);
#endif  

  return;
}

#ifdef FILE_ASCII

  extern void write_properties_ascii(FILE *pfout, struct propiedades_st Prop)
  {
	  type_int   jj;

    fprintf(pfout,"%d %f %f %f ",
         Prop.npart,Prop.pcm[0],Prop.pcm[1],Prop.pcm[2]);
  #ifdef STORE_VELOCITIES  
    fprintf(pfout,"%f %f %f ",
         Prop.vcm[0],Prop.vcm[1],Prop.vcm[2]); 
    fprintf(pfout,"%f %f %f %f %f %f %f ",
         Prop.sig[0],Prop.sig[1],Prop.sig[2],Prop.L[0],Prop.L[1],Prop.L[2],Prop.lambda);
    fprintf(pfout,"%f %f %f %f %f %f %f ",
         Prop.m200,Prop.r200,Prop.v200,Prop.mvir,Prop.rvir,Prop.vvir,Prop.vmax);
    fprintf(pfout,"%f %f ",
         Prop.Ep,Prop.Ec);
  #endif
    fprintf(pfout,"%f %f %f",
         Prop.aa,Prop.bb,Prop.cc);
    for(jj=0;jj<3;jj++)
      fprintf(pfout," %f %f %f",Prop.evec[jj][0],Prop.evec[jj][1],Prop.evec[jj][2]);
  #ifdef STORE_VELOCITIES  
    fprintf(pfout," %f %f %f",
         Prop.aa_vel,Prop.bb_vel,Prop.cc_vel);
    for(jj=0;jj<3;jj++)
      fprintf(pfout," %f %f %f",Prop.evec_vel[jj][0],Prop.evec_vel[jj][1],Prop.evec_vel[jj][2]);
  #endif
    fprintf(pfout,"\n");
    fflush(pfout);
  }

#endif
