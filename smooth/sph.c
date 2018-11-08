#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <math.h>

#include "colores.h"
#include "cosmoparam.h"
#include "variables.h"
#include "sph.h"
#include "grid.h"
#include "colores.h"
#include "bitmask.h"
#include "leesnap.h"

static type_int **sum_fil;
static type_int *id_fil;
static type_int *k_fil;

static void search_cent(const type_real * __restrict__ Pos_cent, type_int *sum_fil, const type_real h2) 
{	
	type_int  id_cent_seg;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  const type_real fac  = (type_real)grid.ngrid/cp.lbox;
  type_real r2, mm, Posprima[3];

  ixc  = (long)(Pos_cent[0]*fac);
  iyc  = (long)(Pos_cent[1]*fac);
  izc  = (long)(Pos_cent[2]*fac);

  ixci = ixc - 1;
  ixcf = ixc + 1;

  iyci = iyc - 1;
  iycf = iyc + 1;

  izci = izc - 1;
  izcf = izc + 1;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= grid.ngrid ) ixcf = grid.ngrid - 1;
  if( iycf >= grid.ngrid ) iycf = grid.ngrid - 1;
  if( izcf >= grid.ngrid ) izcf = grid.ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++)
  {
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= (long)grid.ngrid) ix = ix - (long)grid.ngrid;
    if(ix < 0) ix = ix + grid.ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++)
    {
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= (long)grid.ngrid) iy = iy - (long)grid.ngrid;
      if(iy < 0) iy = iy + (long)grid.ngrid;
      #endif
  
      for( izz = izci ; izz <= izcf ; izz++)
      {
        iz = izz;
        #ifdef PERIODIC
        if(iz >= (long)grid.ngrid) iz = iz - (long)grid.ngrid;
        if(iz < 0) iz = iz + (long)grid.ngrid;
        #endif

        ibox = (ix * (long)grid.ngrid + iy) * (long)grid.ngrid + iz ;
        
        id_cent_seg = grid.icell[ibox];

        while(id_cent_seg != grid.nobj)
        {
          Posprima[0] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]];
          Posprima[1] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]+1];
          Posprima[2] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]+2];

          Posprima[0] -= Pos_cent[0];
          Posprima[1] -= Pos_cent[1];
          Posprima[2] -= Pos_cent[2];

          #ifdef PERIODIC
          if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
          if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
          if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
          if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
          if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
          if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
          #endif

          r2 = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

          if(r2<h2)
          {          
            Posprima[0] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]];
            Posprima[1] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]+1];
            Posprima[2] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]+2];

            Posprima[0] -= Seg[id_fil[id_cent_seg]].Pos_list[0];
            Posprima[1] -= Seg[id_fil[id_cent_seg]].Pos_list[1];
            Posprima[2] -= Seg[id_fil[id_cent_seg]].Pos_list[2];

            #ifdef PERIODIC
            if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
            if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
            if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
            if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
            if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
            if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
            #endif

            mm = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

            if(mm < Seg[id_fil[id_cent_seg]].Rvir[0]) // uso el cuadrado   
            {
              id_cent_seg = grid.ll[id_cent_seg];
              continue;
            }

            Posprima[0] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]];
            Posprima[1] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]+1];
            Posprima[2] = Seg[id_fil[id_cent_seg]].Pos_list[3*k_fil[id_cent_seg]+2];

            Posprima[0] -= Seg[id_fil[id_cent_seg]].Pos_list[3*(Seg[k_fil[id_cent_seg]].size)];  
            Posprima[1] -= Seg[id_fil[id_cent_seg]].Pos_list[3*(Seg[k_fil[id_cent_seg]].size)+1];
            Posprima[2] -= Seg[id_fil[id_cent_seg]].Pos_list[3*(Seg[k_fil[id_cent_seg]].size)+2];

            #ifdef PERIODIC
            if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
            if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
            if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
            if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
            if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
            if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
            #endif

            mm = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

            if(mm < Seg[id_fil[id_cent_seg]].Rvir[1]) // uso el cuadrado   
            {
              id_cent_seg = grid.ll[id_cent_seg];
              continue;
            }

            sum_fil[id_fil[id_cent_seg]]++;

          }// cierra el if
          
          id_cent_seg = grid.ll[id_cent_seg];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;

}

extern void calculo_rho(void)
{
  type_int i, size;
  type_real aux_len, r[3];
  int k;

  aux_len  = cbrt(3.0/(4.0*M_PI*cp.Mpart*(type_real)cp.npart));
  aux_len *= cp.lbox;
  aux_len *= fof[1];
  size = 0;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  fprintf(stdout,"Init SPH\n");
  fflush(stdout);

  for(i=0;i<cp.nseg;i++)
  {
    if(i%10000==0)
    {
      fprintf(stdout,"Seg %u %.4f\n",i,(float)i/(float)cp.nseg);
      fflush(stdout);
    }

    double *esferas, vvv;
    esferas = (double *) malloc(4*(Seg[i].size+2)*sizeof(double));

    Seg[i].Rvir[0] = cbrt(Seg[i].Mass[0])*aux_len;
    Seg[i].Rvir[1] = cbrt(Seg[i].Mass[1])*aux_len;

    for(k=0;k<Seg[i].size;k++)
    {
      r[0] = Seg[i].Pos_list[3*k]  - Seg[i].Pos_list[0];
      r[1] = Seg[i].Pos_list[3*k+1]- Seg[i].Pos_list[1];
      r[2] = Seg[i].Pos_list[3*k+2]- Seg[i].Pos_list[2];

      if(r[0]> 0.5*cp.lbox) Seg[i].Pos_list[3*k]   -= cp.lbox;
      if(r[1]> 0.5*cp.lbox) Seg[i].Pos_list[3*k+1] -= cp.lbox;
      if(r[2]> 0.5*cp.lbox) Seg[i].Pos_list[3*k+2] -= cp.lbox;

      if(r[0]<-0.5*cp.lbox) Seg[i].Pos_list[3*k]   += cp.lbox;
      if(r[1]<-0.5*cp.lbox) Seg[i].Pos_list[3*k+1] += cp.lbox;
      if(r[2]<-0.5*cp.lbox) Seg[i].Pos_list[3*k+2] += cp.lbox;

      esferas[0+4*k] = Seg[i].Pos_list[3*k]  ;
      esferas[1+4*k] = Seg[i].Pos_list[3*k+1];
      esferas[2+4*k] = Seg[i].Pos_list[3*k+2];
      esferas[3+4*k] = R_SPH; 
    }

    k = Seg[i].size;
    esferas[0+4*k] = Seg[i].Pos_list[0];
    esferas[1+4*k] = Seg[i].Pos_list[1];
    esferas[2+4*k] = Seg[i].Pos_list[2];
    esferas[3+4*k] = Seg[i].Rvir[0];

    k = Seg[i].size + 1;
    esferas[0+4*k] = Seg[i].Pos_list[3*(Seg[i].size-1)];
    esferas[1+4*k] = Seg[i].Pos_list[3*(Seg[i].size-1)+1];
    esferas[2+4*k] = Seg[i].Pos_list[3*(Seg[i].size-1)+2];
    esferas[3+4*k] = Seg[i].Rvir[1];
  
    vvv = 0.0;
    k = Seg[i].size + 2;
    arvo_(&k, esferas, &vvv);
    Seg[i].vol = (float)vvv;

    k = 0;
    esferas[0+4*k] = Seg[i].Pos_list[0];
    esferas[1+4*k] = Seg[i].Pos_list[1];
    esferas[2+4*k] = Seg[i].Pos_list[2];
    esferas[3+4*k] = Seg[i].Rvir[0];

    k = 1;
    esferas[0+4*k] = Seg[i].Pos_list[3*(Seg[i].size-1)];
    esferas[1+4*k] = Seg[i].Pos_list[3*(Seg[i].size-1)+1];
    esferas[2+4*k] = Seg[i].Pos_list[3*(Seg[i].size-1)+2];
    esferas[3+4*k] = Seg[i].Rvir[1];

    vvv = 0.0;
    k   = 2; 
    arvo_(&k, esferas, &vvv);
    Seg[i].vol -= (float)vvv;
    Seg[i].vol *= (Seg[i].vol<0.0 ? 0.0f : 1e-9);

    Seg[i].Rvir[0] *= Seg[i].Rvir[0]; // Uso el cuadrado 
    Seg[i].Rvir[1] *= Seg[i].Rvir[1]; //
    size += Seg[i].size;

    free(esferas);
  }

  id_fil  = (type_int *) malloc(size*sizeof(type_int));
  k_fil   = (type_int *) malloc(size*sizeof(type_int));
  sum_fil = (type_int **) calloc(NTHREADS,sizeof(type_int *));
  for(i=0;i<NTHREADS;i++)
    sum_fil[i] = (type_int *) calloc(cp.nseg,sizeof(type_int));

  grid.nobj = size;
  grid.ngrid = (long)(cp.lbox/R_SPH);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }else{
    fprintf(stdout,"Using NGRID = %lu\n",grid.ngrid);
  }

  fflush(stdout);

  grid_init();
  grid_build(id_fil, k_fil);

  read_gadget();

  #pragma omp parallel for schedule(dynamic) num_threads(NTHREADS) \
  default(none) private(i,k) shared(cp,sum_fil,P,Seg,stdout) 
  for(i=0;i<cp.npart;i++)
  {
    if(i%10000000==0)
    {
      fprintf(stdout,"Particles %u %.4f\n",i,(float)i/(float)cp.npart);
      fflush(stdout);
    }

    k = omp_get_thread_num();
    search_cent(P[i].Pos,sum_fil[k],R_SPH*R_SPH);
  }

  #pragma omp parallel for num_threads(NTHREADS) default(none) \
  private(i,k,aux_len,r) shared(cp,Seg,sum_fil,stdout)
  for(i=0;i<cp.nseg;i++)
  {
    Seg[i].Rvir[0] = sqrt(Seg[i].Rvir[0]);
    Seg[i].Rvir[1] = sqrt(Seg[i].Rvir[1]);
    
    Seg[i].mass_part = 0.f;
    for(k=0;k<NTHREADS;k++)
      Seg[i].mass_part += cp.Mpart*(type_real)sum_fil[k][i];

    aux_len = Seg[i].len - Seg[i].Rvir[0] - Seg[i].Rvir[1];
    aux_len *= 1e-3;

    Seg[i].mu  = aux_len<0.0 ? 0.0f : Seg[i].mass_part/aux_len;
    Seg[i].rho = aux_len<0.0 ? 0.0f : Seg[i].mass_part/Seg[i].vol;

    for(k=0;k<Seg[i].size;k++)
    {
      r[0] = Seg[i].Pos_list[3*k+0]+Seg[i].Pos_list[0];
      r[1] = Seg[i].Pos_list[3*k+1]+Seg[i].Pos_list[1];
      r[2] = Seg[i].Pos_list[3*k+2]+Seg[i].Pos_list[2];

      if(r[0]> 0.5*cp.lbox) Seg[i].Pos_list[3*k]   -= cp.lbox;
      if(r[1]> 0.5*cp.lbox) Seg[i].Pos_list[3*k+1] -= cp.lbox;
      if(r[2]> 0.5*cp.lbox) Seg[i].Pos_list[3*k+2] -= cp.lbox;

      if(r[0]<-0.5*cp.lbox) Seg[i].Pos_list[3*k]   += cp.lbox;
      if(r[1]<-0.5*cp.lbox) Seg[i].Pos_list[3*k+1] += cp.lbox;
      if(r[2]<-0.5*cp.lbox) Seg[i].Pos_list[3*k+2] += cp.lbox;
    }
  }

  fprintf(stdout,"End SPH\n");
  fflush(stdout);

  for(k=0;k<NTHREADS;k++)
    free(sum_fil[k]);
  free(sum_fil);
  free(id_fil);
  free(k_fil);
  grid_free();
  free(P);

  return;
}
