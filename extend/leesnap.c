#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

extern void leeheader(void)
{

  FILE *pf;
  char filename[200];
  type_int d1,d2;

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fclose(pf);

  /* Definicion estructura cosmoparam */
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;
  cp.npart     = header.npartTotal[1];
  cp.Mpart     = header.mass[1];
  //cp.Mpart    = 3.143E-4*cp.hparam;  /*A-5*/
  //cp.Mpart    = 3.929E-5*cp.hparam;  /*A-4*/
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

  printf("*********************************** \n");
  printf("*   Parametros de la simulacion   * \n");
  printf("*********************************** \n");
  printf("  Numero de particulas = %u \n", cp.npart);
  printf("  Lado del box = %g \n", cp.lbox);
  printf("  Redshift = %g \n", cp.redshift);
  printf("  Omega Materia = %g \n", cp.omegam);
  printf("  Omega Lambda = %g \n", cp.omegal);
  printf("  Parametro de Hubble = %g \n",cp.hparam);
  printf("  Masa por particula = %g \n",cp.Mpart);
  printf("  Softening = %g\n",cp.soft);
  printf("*********************************** \n");
  printf("*********************************** \n");
}

static void set_directory(const char * prefix, char * name, const type_int NNN, const type_real * fof)
{

  sprintf(name,"./%.2d_%.4d",snap.num,NNN);

  #ifndef ORIGINAL
    sprintf(name,"%s_smooth",name);
  #endif

  #ifdef NEW
    sprintf(name,"%s_new_%s",name,prefix);
  #else
    sprintf(name,"%s_%s",name,prefix);
  #endif

  #ifdef MCRITIC
	  sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);  

	return;
}

static void set_name(const char * prefix, char * name, const type_int NNN, const type_real * fof)
{
  #ifdef ORIGINAL
    sprintf(name,"../");
  #else
    #ifdef NEW_VERSION
      sprintf(name,"../smooth_particles/");
    #else
      sprintf(name,"../smooth/");
    #endif
  #endif

  sprintf(name,"%s%.2d_%.4d",name,snap.num,NNN);

  #ifndef ORIGINAL
    sprintf(name,"%s_smooth",name);
  #endif

  #ifdef NEW
    sprintf(name,"%s_new_%s",name,prefix);
  #else
    sprintf(name,"%s_%s",name,prefix);
  #endif

  #ifdef MCRITIC
	  sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);  

	return;
}

extern void read_segment(type_int NNN, type_real *fof)
{
  char  filename[200];
  type_int   i,j,k,l;
  type_real  r[3];
  type_real  size;
  type_real  *Pos_aux;
  FILE  *pf, *pfprop;

  set_name("segmentos",filename,NNN,fof);
  fprintf(stdout,"%s\n",filename);

  pf = fopen(filename,"rb"); 

  fread(&cp.nseg,sizeof(type_int),1,pf);

  fprintf(stdout,"Segmentos %d\n",cp.nseg);
  fflush(stdout);

  Seg = (struct segmentstd *) malloc(cp.nseg*sizeof(struct segmentstd));

  for(i=0;i<cp.nseg;i++)
  {
    fread(&Seg[i].size,sizeof(type_int),1,pf);
    
    Seg[i].Pos = (type_real *) malloc(3*Seg[i].size*sizeof(type_real));

    for(k=0;k<Seg[i].size;k++)
    {
      #ifdef ORIGINAL
        fread(&j,sizeof(type_int),1,pf);
        memcpy(r[0],Gr[j].Pos[0],3*sizeof(type_real));
      #else
        fread(&r[0], sizeof(type_real), 3, pf);
      #endif

      for(j=0;j<3;j++)
        Seg[i].Pos[3*k+j] = r[j];
    }

  }

  fclose(pf);

  set_name("propiedades",filename,NNN,fof);
  fprintf(stdout,"%s\n",filename);

  pf = fopen(filename,"rb"); 

  fread(&k,sizeof(type_int),1,pf);

  assert(k==cp.nseg);

  fprintf(stdout,"Propiedades Segmentos %d\n",cp.nseg);
  fflush(stdout);

  for(i=0;i<cp.nseg;i++)
  {  
    fread(&Seg[i].flag,sizeof(type_int),1,pf);
    fread(&k,sizeof(type_int),1,pf);
    fread(&Seg[i].Mass[0],sizeof(type_real),2,pf);
    fread(&Seg[i].Vnodos[0],sizeof(type_real),6,pf);
    fread(&Seg[i].razon,sizeof(type_real),1,pf);
    fread(&Seg[i].len,sizeof(type_real),1,pf);
    fread(&Seg[i].elong,sizeof(type_real),1,pf);
    fread(&Seg[i].rms,sizeof(type_real),1,pf);
    #ifdef NEW_VERSION
    fread(&Seg[i].mass_part,sizeof(type_real),1,pf);
    fread(&Seg[i].vol,sizeof(type_real),1,pf);
    fread(&Seg[i].rho,sizeof(type_real),1,pf);
    fread(&Seg[i].mu,sizeof(type_real),1,pf);
    #endif
   
    assert(k==Seg[i].size);
  }

  fclose(pf);

  GREEN("************* Finaliza la lectura *****************\n");
  fflush(stdout);

  GREEN("************* Comienza la escritura *****************\n");
  fflush(stdout);

  set_directory("segmentos_extend",filename,NNN,fof);
  pf = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pf);

  set_directory("propiedades_extend",filename,NNN,fof);
  pfprop = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pfprop);

  for(i=0;i<cp.nseg;i++)
  {
    k = (type_int)(0.5*Seg[i].len/RAUX);
    j = Seg[i].size + 2*k;
    Pos_aux = (type_real *) malloc(3*j*sizeof(type_real));

    fwrite(&j,sizeof(type_int),1,pf);

    #ifdef LAST_VECTOR

      size = 0.0f;
      for(j=0;j<3;j++)
      {
        r[j] = Seg[i].Pos[j+3]-Seg[i].Pos[j];
        #ifdef PERIODIC
        if(r[j]> 0.5*cp.lbox) r[j] -= cp.lbox;
        if(r[j]<-0.5*cp.lbox) r[j] += cp.lbox;
        #endif
        size += (r[j]*r[j]);
      }
      size = sqrt(size);

    #else

      size = 0.0f;
      for(j=0;j<3;j++)
      {
        r[j] = Seg[i].Pos[3*(Seg[i].size-1)+j] - Seg[i].Pos[j];
        #ifdef PERIODIC
        if(r[j]> 0.5*cp.lbox) r[j] -= cp.lbox;
        if(r[j]<-0.5*cp.lbox) r[j] += cp.lbox;
        #endif
        size += (r[j]*r[j]);
      }
      size = sqrt(size);

    #endif

    for(j=0;j<3;j++)
      r[j] /= size;

    for(j=0;j<3;j++)
    {
      Pos_aux[j] = Seg[i].Pos[j] - 0.5*Seg[i].len*r[j];
      #ifdef PERIODIC
      if(Pos_aux[j] > cp.lbox) Pos_aux[j] -= cp.lbox;
      if(Pos_aux[j] <    0.0f) Pos_aux[j] += cp.lbox;
      #endif
    }
  
    size = 0.0f;
    for(l=1;l<k;l++)
    {
      for(j=0;j<3;j++)
      {
        Pos_aux[3*l+j] = Pos_aux[3*(l-1)+j] + RAUX*r[j];
        #ifdef PERIODIC
        if(Pos_aux[3*l+j] >= cp.lbox) Pos_aux[3*l+j] -= cp.lbox;
        if(Pos_aux[3*l+j] <     0.0f) Pos_aux[3*l+j] += cp.lbox;
        #endif
      }
      size += RAUX;
    }
  
    assert(0.5*Seg[i].len-size>0.0f);
    
    ////////////////////////////////////////////
    for(l=0;l<Seg[i].size;l++)
    {
      for(j=0;j<3;j++)
      {
        Pos_aux[3*(l+k)+j] = Seg[i].Pos[3*l+j];
        #ifdef PERIODIC
        if(Pos_aux[3*(l+k)+j] >= cp.lbox) Pos_aux[3*(l+k)+j] -= cp.lbox;
        if(Pos_aux[3*(l+k)+j] <     0.0f) Pos_aux[3*(l+k)+j] += cp.lbox;
        #endif
      }
    }
    ////////////////////////////////////////////
    #ifdef LAST_VECTOR

      size = 0.0f;
      for(j=0;j<3;j++)
      {
        r[j] = Seg[i].Pos[3*(Seg[i].size-1)+j]-Seg[i].Pos[3*(Seg[i].size-2)+j];
        #ifdef PERIODIC
        if(r[j]> 0.5*cp.lbox) r[j] -= cp.lbox;
        if(r[j]<-0.5*cp.lbox) r[j] += cp.lbox;
        #endif
        size += (r[j]*r[j]);
      }
      size = sqrt(size);

      for(j=0;j<3;j++)
        r[j] /= size;

    #endif

    l = Seg[i].size+2*k-1;

    for(j=0;j<3;j++)
    {
      Pos_aux[3*l+j] = Seg[i].Pos[3*(Seg[i].size-1)+j] + 0.5*Seg[i].len*r[j];
      #ifdef PERIODIC
      if(Pos_aux[3*l+j] > cp.lbox) Pos_aux[3*l+j] -= cp.lbox;
      if(Pos_aux[3*l+j] <    0.0f) Pos_aux[3*l+j] += cp.lbox;
      #endif
    }

    size = 0.0f;
    for(l=Seg[i].size+2*k-2;l>=Seg[i].size+k;l--)
    {
      for(j=0;j<3;j++)
      {
        Pos_aux[3*l+j] = Pos_aux[3*(l+1)+j] - RAUX*r[j];
        #ifdef PERIODIC
        if(Pos_aux[3*l+j] >= cp.lbox) Pos_aux[3*l+j] -= cp.lbox;
        if(Pos_aux[3*l+j] <     0.0f) Pos_aux[3*l+j] += cp.lbox;
        #endif
      }
      size += RAUX;
    }

    assert(0.5*Seg[i].len-size>0.0f);

    l = Seg[i].size+2*k;
    fwrite(Pos_aux,sizeof(type_real),3*l,pf);

    fwrite(&Seg[i].flag,sizeof(type_int),1,pfprop);
    fwrite(&l,sizeof(type_int),1,pfprop);
    fwrite(Seg[i].Mass,sizeof(type_real),2,pfprop);
    fwrite(Seg[i].Vnodos,sizeof(type_real),6,pfprop);
    fwrite(&Seg[i].razon,sizeof(type_real),1,pfprop);
    fwrite(&Seg[i].len,sizeof(type_real),1,pfprop);
    fwrite(&Seg[i].elong,sizeof(type_real),1,pfprop);
    fwrite(&Seg[i].rms,sizeof(type_real),1,pfprop);
    #ifdef NEW_VERSION
    fwrite(&Seg[i].mass_part,sizeof(type_real),1,pfprop);
    fwrite(&Seg[i].vol,sizeof(type_real),1,pfprop);
    fwrite(&Seg[i].rho,sizeof(type_real),1,pfprop);
    fwrite(&Seg[i].mu,sizeof(type_real),1,pfprop);
    #endif

    free(Pos_aux);

  }

  fclose(pf);
  fclose(pfprop);

  GREEN("************* Finaliza la escritura *****************\n");
  fflush(stdout);

  #ifdef ORIGINAL
  free(Gr);
  #endif

  return;
}

#ifdef ORIGINAL

  extern void read_grup_fof(type_real *fof)
  {
    char  filename[200];
    type_int   i;
    FILE  *pfin;
   
    #ifdef MCRITIC
      sprintf(filename,"../%.2d_%.2f_centros_cut_%.2f.bin",snap.num,fof[1],m_critica);
    #else
      sprintf(filename,"../../%.2d_%.2f_centros.bin",snap.num,fof[1]);
    #endif
  
    pfin = fopen(filename,"rb"); 
  
    fread(&cp.ngrup,sizeof(type_int),1,pfin);
  
    fprintf(stdout,"Grupos %d\n",cp.ngrup);
    fflush(stdout);
  
    Gr = (struct grup_data *) malloc(cp.ngrup*sizeof(struct grup_data));
  
    for(i=0;i<cp.ngrup;i++)
    {
      fread(&Gr[i].save,sizeof(type_int),1,pfin);
      fread(&Gr[i].id,sizeof(type_int),1,pfin);
      fread(&Gr[i].Pos[0],sizeof(float),1,pfin);
      fread(&Gr[i].Pos[1],sizeof(float),1,pfin);
      fread(&Gr[i].Pos[2],sizeof(float),1,pfin);
      fread(&Gr[i].NumPart,sizeof(type_int),1,pfin);
    }
  
    fclose(pfin);
  
    return;
  
  }

#endif
