#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "create_random.h"
#include "colores.h"
#include "leesnap.h"

void create_rand(int NNN, type_real *fof)
{
  char filename[200];
  int i, j, k, idim, Tid, *c;
  int totfil = 0;
  type_real  vdir[3];
  type_real  xc,yc,zc;
  FILE **pfrand;

  #ifdef MPC 
  type_real POSFACTOR = 1000.;
  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(filename,"Reescala lbox %g Mpc to %g Kpc\n",cp.lbox/POSFACTOR,cp.lbox);GREEN(filename);
  GREEN("**********************************\n");
  #endif

  srand48(80); // guarda la semilla
  c = (int *) malloc(NTHREADS*sizeof(int));
  pfrand = (FILE **) malloc(NTHREADS*sizeof(FILE));

  for(i=0;i<NTHREADS;i++)
  {
    c[i] = 0;

    #ifdef MCRITIC
      sprintf(filename,"%.2d_%.4d_random_fil_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NNN,m_critica,fof[0],fof[1],i);
    #else
      sprintf(filename,"%.2d_%.4d_random_fil_%.2f_%.2f.%.2d.bin",snap.num,NNN,fof[0],fof[1],i);
    #endif
    pfrand[i]=fopen(filename,"w");
    fwrite(&c[i],sizeof(int),1,pfrand[i]); // Nfil 
    fwrite(&c[i],sizeof(int),1,pfrand[i]); // Nrand
  }

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(static) default(none) private(i,j,k,idim,Tid, \
  vdir,xc,yc,zc) shared(pfrand,Gr,Seg,cp,c, nrand,stdout) \
  reduction(+:totfil)
  for(i=0;i<cp.nseg;i++)
  {
    Tid = omp_get_thread_num();

  	type_real xang[nrand];
  	type_real yang[nrand];
  	type_real zang[nrand];
    type_real *raux;
    #ifdef RECTA
      type_real **vdir_rand;
    #else
      type_real ***vdir_rand;
    #endif

    ////////////////////////////////////////////////////////////////////
    raux = (type_real *) malloc((Seg[i].size-1)*sizeof(type_real));
    #ifdef RECTA
      vdir_rand = (type_real **) malloc(nrand*sizeof(type_real *));
    #else
      vdir_rand = (type_real ***) malloc(nrand*sizeof(type_real **));
    #endif
    ////////////////////////////////////////////////////////////////////

    for(j=0;j<nrand;j++)
    {
   	  xang[j]  = drand48()*2.0*M_PI;
  	  yang[j]  = drand48()*2.0*M_PI;
  	  zang[j]  = drand48()*2.0*M_PI;

      #ifdef RECTA
        vdir_rand[j] = (type_real *) malloc(3*sizeof(type_real));
      #else
        vdir_rand[j] = (type_real **) malloc((Seg[i].size-1)*sizeof(type_real *));
      #endif
    } 

    for(k=1;k<Seg[i].size;k++)
    {
      raux[k-1] = 0.0;
      for(idim=0;idim<3;idim++)
      {
        vdir[idim] = Gr[Seg[i].list[k]].Pos[idim]-Gr[Seg[i].list[k-1]].Pos[idim];    
        #ifdef PERIODIC
        if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
        if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
        #endif
        raux[k-1] += vdir[idim]*vdir[idim];
      }

      raux[k-1] = sqrt(raux[k-1]);

      for(idim=0;idim<3;idim++)
        vdir[idim] /= raux[k-1];

      #ifdef RECTA
      if(k!=1) continue;
      #endif
 
      for(j=0;j<nrand;j++) 
      {
        #ifdef RECTA
          memcpy(vdir_rand[j],vdir,3*sizeof(type_real));
          rotacion(vdir_rand[j],xang[j],yang[j],zang[j]);
        #else
          vdir_rand[j][k-1] = (type_real *) calloc(3,sizeof(type_real));
          memcpy(vdir_rand[j][k-1],vdir,3*sizeof(type_real));
          rotacion(vdir_rand[j][k-1],xang[j],yang[j],zang[j]);
        #endif
      }
    }

    fwrite(&i,sizeof(int),1,pfrand[Tid]);
    fwrite(&Seg[i].flag,sizeof(int),1,pfrand[Tid]);
    fwrite(&Seg[i].size,sizeof(int),1,pfrand[Tid]);
    fwrite(&Seg[i].len,sizeof(type_real),1,pfrand[Tid]);
    fwrite(&Gr[Seg[i].list[0]].NumPart,sizeof(int),1,pfrand[Tid]);
    fwrite(&Gr[Seg[i].list[Seg[i].size-1]].NumPart,sizeof(int),1,pfrand[Tid]);
    
    for(j=0;j<nrand;j++) 
    {
      fwrite(&xang[j],sizeof(type_real),1,pfrand[Tid]);
      fwrite(&yang[j],sizeof(type_real),1,pfrand[Tid]);
      fwrite(&zang[j],sizeof(type_real),1,pfrand[Tid]);

      xc = Gr[Seg[i].list[0]].Pos[0];
      yc = Gr[Seg[i].list[0]].Pos[1];
      zc = Gr[Seg[i].list[0]].Pos[2];

      fwrite(&xc,sizeof(type_real),1,pfrand[Tid]);
      fwrite(&yc,sizeof(type_real),1,pfrand[Tid]);
      fwrite(&zc,sizeof(type_real),1,pfrand[Tid]);

      for(k=1;k<Seg[i].size;k++)
      {
        #ifdef RECTA
          xc += raux[k-1]*vdir_rand[j][0];
          yc += raux[k-1]*vdir_rand[j][1];
          zc += raux[k-1]*vdir_rand[j][2];   
        #else
          xc += raux[k-1]*vdir_rand[j][k-1][0];
          yc += raux[k-1]*vdir_rand[j][k-1][1];
          zc += raux[k-1]*vdir_rand[j][k-1][2];   
        #endif

        #ifdef PERIODIC
        xc = xc >= cp.lbox ? xc-cp.lbox : xc;
        xc = xc <      0.0 ? xc+cp.lbox : xc;
      
        yc = yc >= cp.lbox ? yc-cp.lbox : yc;
        yc = yc <      0.0 ? yc+cp.lbox : yc;
 
        zc = zc >= cp.lbox ? zc-cp.lbox : zc;
        zc = zc <      0.0 ? zc+cp.lbox : zc;
        #endif

        fwrite(&xc,sizeof(type_real),1,pfrand[Tid]);
        fwrite(&yc,sizeof(type_real),1,pfrand[Tid]);
        fwrite(&zc,sizeof(type_real),1,pfrand[Tid]);

        #ifndef RECTA
          free(vdir_rand[j][k-1]);
        #endif

      }
      free(vdir_rand[j]);
    }
    
    free(vdir_rand);

    free(raux);
    c[Tid]++;
    totfil++;

  } //FINALIZA EL PARALELO

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,nrand*c[i]);  
    rewind(pfrand[i]);
    fwrite(&c[i],sizeof(int),1,pfrand[i]);
    fwrite(&nrand,sizeof(int),1,pfrand[i]);
    fclose(pfrand[i]);
  }

  free(c);

  BLUE("******************************\n");
  fprintf(stdout,"Calcula %d Filamentos\n",totfil);  
  fprintf(stdout,"Calcula %d Filamentos Random\n",nrand*totfil);  
  BLUE("******************************\n");

  return;
}

void rotacion(type_real dd[3], type_real xphi, type_real yphi, type_real zphi)
{
  type_real aux_vec[3], vec[3];

// rota alrededor de Z
// un angulo xphi
  vec[0] =  dd[0]*cos(zphi)+dd[1]*sin(zphi);
  vec[1] = -dd[0]*sin(zphi)+dd[1]*cos(zphi);
  vec[2] =  dd[2];

// rota alrededor de Y
// un angulo yphi
  aux_vec[0] =  vec[0]*cos(yphi)-vec[2]*sin(yphi);
  aux_vec[1] =  vec[1];
  aux_vec[2] =  vec[0]*sin(yphi)+vec[2]*cos(yphi);

// rota alrededor de X
// un angulo xphi 
  dd[0] = aux_vec[0];
  dd[1] = aux_vec[1]*cos(xphi)+aux_vec[2]*sin(xphi);
  dd[2] =-aux_vec[1]*sin(xphi)+aux_vec[2]*cos(xphi);
 
  return;
}
