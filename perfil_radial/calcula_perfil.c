#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "list.h"
#include "colores.h"
#include "leesnap.h"
#include "calcula_perfil.h"

void propiedades_halos(type_real *fof)
{
  char filename[200];
  int i, j, k, idim, Tid, *c;
  int tothalos = 0;
  int binsper = nbins + 1;
  type_real *rsph2;                 // En Kpc                    
  type_real *volinv;                // En Kpc
  type_real r;                      // En Kpc
  type_real delta;                 
  type_real RMAX = 20000., RMIN = RMAX/100.;
  FILE **pfdens;
  #ifdef CALCULA_MEDIA
  int *npart;
  FILE **pfvel;
  #endif
  #ifdef MPC 
  type_real POSFACTOR = 1000.;

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(filename,"Reescala lbox %g Mpc to %g Kpc\n",cp.lbox/POSFACTOR,cp.lbox);GREEN(filename);
  GREEN("**********************************\n");
  #endif

  fprintf(stdout,"Using Nbins = %d\n",nbins);
  fprintf(stdout,"RMAX %f Mpc\n",RMAX/1000.);
  fprintf(stdout,"RMIN %f Mpc\n",RMIN/1000.);

  //////////////////////////////////////////////////
  fprintf(stdout,"Build grid\n");

  grid.nobj = cp.npart;
  grid.step = 1000.;
  //grid.step = RMAX;
  grid.ngrid = (int)(cp.lbox/grid.step);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
    grid.step = cp.lbox/(type_real)grid.ngrid;
  }

  grid_init();
  grid_build();
  //////////////////////////////////////////////////

  c = (int *) malloc(NTHREADS*sizeof(int));
  pfdens = (FILE **) malloc(NTHREADS*sizeof(FILE));
  #ifdef CALCULA_MEDIA
  pfvel = (FILE **) malloc(NTHREADS*sizeof(FILE));
  #endif

  for(i=0;i<NTHREADS;i++)
  {
    c[i] = 0;

    #ifdef MCRITIC
    sprintf(filename,"%.2d_halo_densidad_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,m_critica,fof[0],fof[1],i);
    #else
    sprintf(filename,"%.2d_halo_densidad_%.2f_%.2f.%.2d.bin",snap.num,fof[0],fof[1],i);
    #endif
    pfdens[i]=fopen(filename,"w");
    fwrite(&c[i],sizeof(int),1,pfdens[i]);        

    #ifdef CALCULA_MEDIA
      #ifdef MCRITIC
        sprintf(filename,"../%.2d_halo_vmedia_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,m_critica,fof[0],fof[1],i);
      #else
        sprintf(filename,"../%.2d_halo_vmedia_%.2f_%.2f.%.2d.bin",snap.num,fof[0],fof[1],i);
      #endif
    pfvel[i]=fopen(filename,"w");
    fwrite(&c[i],sizeof(int),1,pfvel[i]);        
    fwrite(&nbins,sizeof(int),1,pfvel[i]);
    #endif
  }


  #ifdef CALCULA_MEDIA

  BLUE("**** CALCULA VEL MEDIAS ******\n");

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) private(i,j,k, \
  rsph2,npart,Tid) shared(pfdens,Gr,cp,RMAX,RMIN, \
  c,binsper,nbins,stdout) 
  for(i=0;i<cp.ngrup;i++)
  {
    Tid = omp_get_thread_num();
    npart = (int *) calloc(nbins,sizeof(int));
    Gr[i].Vmedia = (type_real *) calloc(3*nbins,sizeof(type_real));
    rsph2  = (type_real *) malloc(binsper*sizeof(type_real));

    #ifdef LOG_BINS
      logspace(rsph2,RMAX,RMIN,binsper);
    #else
      linspace(rsph2,RMAX,RMIN,binsper);
    #endif

    calc_media(Gr[i].Pos[0],Gr[i].Pos[1],Gr[i].Pos[2],npart,Gr[i].Vmedia,rsph2,nbins);

    for(j=0;j<nbins;j++)
    {
      if(npart[j]==0)
        continue;

      for(k=0;k<3;k++)
        Gr[i].Vmedia[j+nbins*k] /= (type_real)npart[j];
    }

    free(npart); 
    free(rsph2);

    c[Tid]++;

  } //FINALIZA EL PARALELO

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfvel[i]);
    fwrite(&c[i],sizeof(int),1,pfvel[i]);
    fclose(pfvel[i]);
    c[i] = 0;
  }

  #endif
 
  BLUE("********** CALCULA ***********\n"); 
  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) private(i,j,k,idim, \
  rsph2,Tid,volinv,delta,r) shared(pfdens,Gr,cp,RMAX,RMIN, \
  c,binsper,nbins,stdout) reduction(+:tothalos)
  for(i=0;i<cp.ngrup;i++)
  {
    int        *vsize;
    type_real  aux_mean, aux_rms;
    type_real **vmean, **vquad;
    type_real *vmean_rad, *vquad_rad;
    #ifdef HUECOS
    type_real tmpvol;
    #endif

    rsph2  = (type_real *) malloc(binsper*sizeof(type_real));
    volinv = (type_real *) malloc(binsper*sizeof(type_real));

    #ifdef LOG_BINS
      logspace(rsph2,RMAX,RMIN,binsper);
    #else
      linspace(rsph2,RMAX,RMIN,binsper);
    #endif

    for(k=binsper-1;k>=0;k--)
    {
      #ifdef HUECOS
        volinv[k] = 4.0*M_PI*sqrt(rsph2[k])*rsph2[k]/3.0;

        if(k!=binsper-1) 
          volinv[k] -= tmpvol;        

        tmpvol = 4.0*M_PI*sqrt(rsph2[k])*rsph2[k]/3.0;
      #endif  
        volinv[k] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[k])/(type_real)cp.npart;
    }

    Tid = omp_get_thread_num();

    vsize     = (int *) malloc(nbins*sizeof(int));
    vmean_rad = (type_real *) malloc(nbins*sizeof(type_real));
    vquad_rad = (type_real *) malloc(nbins*sizeof(type_real));

    vmean     = (type_real **) malloc(nbins*sizeof(type_real *));
    vquad     = (type_real **) malloc(nbins*sizeof(type_real *));
    for(j=0;j<nbins;j++)
    {
      vmean[j] = (type_real *) malloc(3*sizeof(type_real));
      vquad[j] = (type_real *) malloc(3*sizeof(type_real));
    }

    calc_part(Gr[i].Pos[0],Gr[i].Pos[1],Gr[i].Pos[2],i,vsize,vmean,vquad,vmean_rad,vquad_rad,rsph2,nbins);

    fwrite(&i,sizeof(int),1,pfdens[Tid]);                   // Id
    fwrite(&Gr[i].NumPart,sizeof(int),1,pfdens[Tid]);       // escribe el numero de particulas 
    r = Gr[i].NumPart*cp.Mpart;
    fwrite(&r,sizeof(type_real),1,pfdens[Tid]);             // escribe la masa

    fwrite(&nbins,sizeof(int),1,pfdens[Tid]);                // la cantidad de pelotas que hizo

    for(j=0;j<nbins;j++)
    {
      //r  = sqrt(rsph2[j]/rsph2[0]);
      //r += sqrt(rsph2[j+1]/rsph2[0]);
      r  = sqrt(rsph2[j]);
      r += sqrt(rsph2[j+1]);
      r *= 0.5;

      fwrite(&r,sizeof(type_real),1,pfdens[Tid]);       // rbin_centre

      delta = (type_real)vsize[j]*volinv[j]-1.0;

      fwrite(&vsize[j],sizeof(int),1,pfdens[Tid]);       // npart
      fwrite(&delta,sizeof(type_real),1,pfdens[Tid]);    // pho - overdense

      ////////////////////////////////////////////////////////////////////////////////////

      if(vsize[j]!=0)
      {
        for(idim=0;idim<3;idim++)
        {
          vmean[j][idim] /= (type_real)vsize[j];
          vquad[j][idim] /= (type_real)vsize[j];
        }
        vmean_rad[j] /= (type_real)vsize[j];
        vquad_rad[j] /= (type_real)vsize[j];
      }

      aux_mean = aux_rms = 0.0;

      for(idim=0;idim<3;idim++)
      {
        aux_mean  += vmean[j][idim];
        aux_rms   += vquad[j][idim] - vmean[j][idim]*vmean[j][idim];
      }

      aux_mean /= 3.0;
      if(vsize[j]>10)
        aux_rms  = sqrt(fabs(aux_rms))/3.0;
      else
        aux_rms  = 0.0;

      fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean
      fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms

      aux_mean = vmean_rad[j];
      if(vsize[j]>10)
        aux_rms  = sqrt(fabs(vquad_rad[j] - vmean_rad[j]*vmean_rad[j]));
      else
        aux_rms  = 0.0;

      fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean par
      fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms  par

    }

    for(j=0;j<nbins;j++)
    {
      free(vmean[j]);
      free(vquad[j]);
    }

    free(vmean);
    free(vquad); 
    free(vmean_rad);
    free(vquad_rad); 
    free(rsph2);
    free(volinv);

    c[Tid]++;
    tothalos++;

  } //FINALIZA EL PARALELO

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfdens[i]);
    fwrite(&c[i],sizeof(int),1,pfdens[i]);
    fclose(pfdens[i]);
  }

  free(c);
  grid_free();

  BLUE("******************************\n");
  fprintf(stdout,"Calcula %d Halos\n",tothalos);  
  BLUE("******************************\n");

  return;
}

void calc_part(type_real xc, type_real yc, type_real zc, int icent, int *npart, type_real **mean, type_real **quad, \
               type_real *mean_rad, type_real *quad_rad, type_real *rsph2, int nsph)
{
  int i,j,k;
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  type_real lbox,fac,lbox2;
  type_real dis;
  type_real aux_dot, aux_mod;
  int ngrid, ilim;
  type_real Posprima[3],Vprima[3];
  #ifdef HUECOS
  int bin;
  #endif

  ////////////////////////////////////////////////////
  memset(npart,0,nsph*sizeof(int));
  memset(mean_rad,0.0,nsph*sizeof(type_real));
  memset(quad_rad,0.0,nsph*sizeof(type_real));
  for(j=0;j<nsph;j++)
  {
    memset(mean[j],0.0,3*sizeof(type_real));
    memset(quad[j],0.0,3*sizeof(type_real));
  }
  ////////////////////////////////////////////////////

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  ilim  = (int)(sqrt(2.*rsph2[0])*fac)+1;

  ixc  = (int)(xc*fac);
  ixci = ixc - ilim;
  ixcf = ixc + ilim;
  iyc  = (int)(yc*fac);
  iyci = iyc - ilim;
  iycf = iyc + ilim;
  izc  = (int)(zc*fac);
  izci = izc - ilim;
  izcf = izc + ilim;

  //if(icent==10)
  //{
  //  printf("%d\n",icent);
  //  printf("%d %d %d\n",ixc,ixci,ixcf);
  //  printf("%d %d %d\n",iyc,iyci,iycf);
  //  printf("%d %d %d\n",izc,izci,izcf);
  //}

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= ngrid ) ixcf = ngrid - 1;
  if( iycf >= ngrid ) iycf = ngrid - 1;
  if( izcf >= ngrid ) izcf = ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= ngrid) ix = ix - ngrid;
    if(ix < 0) ix = ix + ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= ngrid) iy = iy - ngrid;
      if(iy < 0) iy = iy + ngrid;
      #endif
  
      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= ngrid) iz = iz - ngrid;
        if(iz < 0) iz = iz + ngrid;
        #endif

        ibox = (ix * ngrid + iy) * ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {

          Posprima[0] = P[i].Pos[0] - xc;
          Posprima[1] = P[i].Pos[1] - yc;
          Posprima[2] = P[i].Pos[2] - zc;

          #ifdef PERIODIC
          if(Posprima[0] >  lbox2) Posprima[0] = Posprima[0] - lbox;
          if(Posprima[1] >  lbox2) Posprima[1] = Posprima[1] - lbox;
          if(Posprima[2] >  lbox2) Posprima[2] = Posprima[2] - lbox;
          if(Posprima[0] < -lbox2) Posprima[0] = Posprima[0] + lbox;
          if(Posprima[1] < -lbox2) Posprima[1] = Posprima[1] + lbox;
          if(Posprima[2] < -lbox2) Posprima[2] = Posprima[2] + lbox;
          #endif

          dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

          #ifdef HUECOS

          //if(dis<rsph2[0] && dis>rsph2[nsph])
          if(dis<rsph2[0])
          {

            bin = 0;
            for(j=1;j<nsph;j++)
            {
              if(dis<rsph2[j])
               bin = j;
              else
                break;
            }          

            aux_mod = sqrt(dis);
            
            aux_dot = 0.0;
            for(k=0;k<3;k++)
            {
              Vprima[k] = P[i].Vel[k]-Gr[icent].Vmedia[bin+nsph*k];
              //Vprima[k] = P[i].Vel[k];

              mean[bin][k] += Vprima[k];
              quad[bin][k] += (Vprima[k]*Vprima[k]);
    
              Posprima[k] /= aux_mod;

              aux_dot += (Vprima[k]*Posprima[k]);
            }
  
            mean_rad[bin] += aux_dot;
            quad_rad[bin] += (aux_dot*aux_dot);
   
            npart[bin]++;
          }// cierra el if

          #endif

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

#ifdef CALCULA_MEDIA

void calc_media(type_real xc, type_real yc, type_real zc, int *numpart, type_real *vel_media, type_real *rsph2, int binsper)
{
  int i,j,k;
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  type_real Posprima[3];
  type_real lbox,fac,lbox2;
  type_real dis;
  int ngrid, ilim;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  ilim  = (int)(sqrt(2.*rsph2[0])*fac)+1;

  ixc  = (int)(xc*fac);
  ixci = ixc - ilim;
  ixcf = ixc + ilim;
  iyc  = (int)(yc*fac);
  iyci = iyc - ilim;
  iycf = iyc + ilim;
  izc  = (int)(zc*fac);
  izci = izc - ilim;
  izcf = izc + ilim;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= ngrid ) ixcf = ngrid - 1;
  if( iycf >= ngrid ) iycf = ngrid - 1;
  if( izcf >= ngrid ) izcf = ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= ngrid) ix = ix - ngrid;
    if(ix < 0) ix = ix + ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= ngrid) iy = iy - ngrid;
      if(iy < 0) iy = iy + ngrid;
      #endif
  
      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= ngrid) iz = iz - ngrid;
        if(iz < 0) iz = iz + ngrid;
        #endif

        ibox = (ix * ngrid + iy) * ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {
          Posprima[0] = P[i].Pos[0] - xc;
          Posprima[1] = P[i].Pos[1] - yc;
          Posprima[2] = P[i].Pos[2] - zc;

          #ifdef PERIODIC
          if(Posprima[0] >  lbox2) Posprima[0] = Posprima[0] - lbox;
          if(Posprima[1] >  lbox2) Posprima[1] = Posprima[1] - lbox;
          if(Posprima[2] >  lbox2) Posprima[2] = Posprima[2] - lbox;
          if(Posprima[0] < -lbox2) Posprima[0] = Posprima[0] + lbox;
          if(Posprima[1] < -lbox2) Posprima[1] = Posprima[1] + lbox;
          if(Posprima[2] < -lbox2) Posprima[2] = Posprima[2] + lbox;
          #endif

          dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

          //if(dis<rsph2[0] && dis>rsph2[binsper])
          if(dis<rsph2[0])
          {

           for(k=0;k<3;k++)
             vel_media[binsper*k] += P[i].Vel[k];

           numpart[0]++;

           for(j=1;j<binsper;j++)
           {
              if(dis<rsph2[j])
              {
                for(k=0;k<3;k++)
                  vel_media[j+binsper*k] += P[i].Vel[k];

                numpart[j]++;

              }else{

                break;

              }
            }
          } // cierra dis

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

#endif

int point_inside(type_real dot, type_real rlong_2)
{
  if(fabs(dot)>rlong_2){

      return 0;

  }else{

    return 1;

  }
}

#ifdef LOG_BINS

  void logspace(type_real *rsph2, type_real max, type_real min, int bins) 
  {

    int i;
    type_real end   = log10(max);
    type_real start = log10(min);
    type_real delta = (end - start) / (bins-1);

    rsph2[0] = pow(10.0,end);   
    rsph2[0] *= rsph2[0];

    for(i=1;i<bins; i++)
    {
      rsph2[i] = pow(10.0,(start + delta * (bins - 1 - i)));
      rsph2[i] *= rsph2[i];
    }  

    return;
  }

#else

  void linspace(type_real *rsph2, type_real max, type_real min, int bins) 
  {
  
    int i;
    type_real end   = max;
    type_real start = min;
    type_real delta = (end - start) / (bins-1);
  
    rsph2[0] = end;   
    rsph2[0] *= rsph2[0];
  
    for(i=1;i<bins; i++)
    {
      rsph2[i] = start + delta * (bins - 1 - i);
      rsph2[i] *= rsph2[i];
    }  
  
    return;
  }

#endif
