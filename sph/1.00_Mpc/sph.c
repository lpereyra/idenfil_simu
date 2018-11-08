#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include "colores.h"
#include "cosmoparam.h"
#include "variables.h"
#include "sph.h"
#include "grid.h"
#include "colores.h"

#define  K_NEIGH_SIZE    32
#define  SET_MAX_SIZE    4*K_NEIGH_SIZE
#define  KERNEL_COEFF_1  2.546479089470      // Coefficients for SPH spline kernel and its derivative  
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786      // Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 

static type_real pmin[3];  // Min coordenadas
static type_real pmax[3];  // Max coordenadas
static type_real r00[3], r01[3];

extern void change_positions()
{
  type_int i;

  for(i=0;i<3;i++)
  {
    pmin[i] = 1.E26; 
    pmax[i] = -1.E26;
  }

  for(i=0;i<cp.npart;i++)
  {
    if(P[i].Pos[0] > pmax[0]) pmax[0] = P[i].Pos[0];
    if(P[i].Pos[0] < pmin[0]) pmin[0] = P[i].Pos[0];
    if(P[i].Pos[1] > pmax[1]) pmax[1] = P[i].Pos[1];
    if(P[i].Pos[1] < pmin[1]) pmin[1] = P[i].Pos[1];
    if(P[i].Pos[2] > pmax[2]) pmax[2] = P[i].Pos[2];
    if(P[i].Pos[2] < pmin[2]) pmin[2] = P[i].Pos[2];
    #ifndef SPH
    P[i].flag = 0;
    #endif
  }

  //fprintf(stdout,"xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  //fprintf(stdout,"ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  //fprintf(stdout,"zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);
 
  for(i=0;i<cp.npart;i++)
  {
    P[i].Pos[0] -= pmin[0];
    P[i].Pos[1] -= pmin[1];
    P[i].Pos[2] -= pmin[2];
  }

  cp.lbox = 0.0f;
  for(i = 0; i < 3; i++)
    if(cp.lbox < (pmax[i] - pmin[i])) cp.lbox = (pmax[i] - pmin[i]);

  cp.lbox *= 1.001;

  //sprintf(message,"lbox %f\n",cp.lbox);
  //GREEN(message);

  return;

}

#ifdef SPH
static void search_cent(const type_int ic, const type_real h2) 
#else
static void search_cent(const type_int iseg, const type_real * __restrict__ Pos_cent, const type_real h2) 
#endif
{	
	type_int  i,idim;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  const type_real fac  = (type_real)grid.ngrid/cp.lbox;
  type_real r2, m00, m01, Posprima[3];
  #ifdef SPH

    const type_real hinv = 1.0 / HSPH;
    const type_real hinv3 = hinv * hinv * hinv;
    type_real r, u, wk, mass_j;
    P[ic].rho = P[ic].mass_part = 0.0f;

    ixc  = (long)(P[ic].Pos[0]*fac);
    iyc  = (long)(P[ic].Pos[1]*fac);
    izc  = (long)(P[ic].Pos[2]*fac);

  #else

    ixc  = (long)(Pos_cent[0]*fac);
    iyc  = (long)(Pos_cent[1]*fac);
    izc  = (long)(Pos_cent[2]*fac);

  #endif

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

        i = grid.llirst[ibox];

        while(i != grid.nobj)
        {
          #ifdef SPH
            Posprima[0] = P[i].Pos[0] - P[ic].Pos[0];
            Posprima[1] = P[i].Pos[1] - P[ic].Pos[1];
            Posprima[2] = P[i].Pos[2] - P[ic].Pos[2];
          #else
            Posprima[0] = P[i].Pos[0] - Pos_cent[0];
            Posprima[1] = P[i].Pos[1] - Pos_cent[1];
            Posprima[2] = P[i].Pos[2] - Pos_cent[2];
          #endif

          //#ifdef PERIODIC
          if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
          if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
          if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
          if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
          if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
          if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
          //#endif

          r2 = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

          if(r2<h2)
          {          
            #ifdef SPH

              r = sqrt(r2);

              u = r * hinv;

              if(u < 0.5)
              {
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              }else{  
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
              }
          
              mass_j = 1.0f;
            
              P[ic].rho += mass_j * wk;
          
              P[ic].mass_part += NORM_COEFF * wk / hinv3;

            #else

            m00 = m01 = 0.0f;
            for(idim=0;idim<3;idim++)
            {
              Posprima[idim] = P[i].Pos[idim] - r00[idim];
              if(Posprima[idim] >= 0.5f*cp.lbox) Posprima[idim] -= cp.lbox;
              if(Posprima[idim] < -0.5f*cp.lbox) Posprima[idim] += cp.lbox;
              m00 += Posprima[idim]*Posprima[idim];

              Posprima[idim] =  P[i].Pos[idim] - r01[idim];
              if(Posprima[idim] >= 0.5f*cp.lbox) Posprima[idim] -= cp.lbox;
              if(Posprima[idim] < -0.5f*cp.lbox) Posprima[idim] += cp.lbox;
              m01 += Posprima[idim]*Posprima[idim];
            }
            m00 = sqrt(m00);
            m01 = sqrt(m01);

            if(m00>=Seg[iseg].Rvir[0] && m01>=Seg[iseg].Rvir[1])              
              P[i].flag = 1;

            #endif
          }// cierra el if

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;

}
 
extern void init_sph(const type_int i)
{
  int j;
  type_int k;
  type_int vec = 0;
  type_real r[3], fac_len;
  const type_real h2  = HSPH*HSPH;
  double *x_pos, *y_pos, *z_pos, *r_sph;
  double vvv = 0.0f;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  for(j=0;j<3;j++)
  {
    r00[j] = Seg[i].Pos[j]-pmin[j];
    if(r00[j]> 0.5*cp.simulbox) r00[j] -= cp.simulbox;
    if(r00[j]<-0.5*cp.simulbox) r00[j] += cp.simulbox;

    r01[j] = Seg[i].Pos[3*(Seg[i].size-1)+j]-pmin[j];
    if(r01[j]> 0.5*cp.simulbox) r01[j] -= cp.simulbox;
    if(r01[j]<-0.5*cp.simulbox) r01[j] += cp.simulbox;
  }

  j = Seg[i].size+2;
  x_pos = (double *) malloc(j*sizeof(double));
  y_pos = (double *) malloc(j*sizeof(double));
  z_pos = (double *) malloc(j*sizeof(double));
  r_sph = (double *) malloc(j*sizeof(double));

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) private(j,k,r) \
  shared(cp,Seg,pmin,x_pos,y_pos,z_pos,r_sph,stdout)
  for(k=0;k<Seg[i].size;k++)
  {
    
    for(j=0;j<3;j++)
    {
      r[j] = Seg[i].Pos[3*k+j]-pmin[j];
      if(r[j]> 0.5*cp.simulbox) r[j] -= cp.simulbox;
      if(r[j]<-0.5*cp.simulbox) r[j] += cp.simulbox;
    }

    #ifdef SPH
      search_cent(i,h2);
    #else
      search_cent(i,r,h2);
    #endif

    x_pos[k] = Seg[i].Pos[3*k];
    y_pos[k] = Seg[i].Pos[3*k+1];
    z_pos[k] = Seg[i].Pos[3*k+2];
    r_sph[k] = HSPH; 

  }

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) private(k) \
  shared(P,cp) reduction(+:vec)
  for(k=0;k<cp.npart;k++)
  {
    if(P[k].flag == 1)
      vec++;
  }

  k = Seg[i].size;

  x_pos[k] = Seg[i].Pos[0];
  y_pos[k] = Seg[i].Pos[1];
  z_pos[k] = Seg[i].Pos[2];
  r_sph[k] = Seg[i].Rvir[0];
  
  k = Seg[i].size + 1;

  x_pos[k] = Seg[i].Pos[3*(Seg[i].size-1)];
  y_pos[k] = Seg[i].Pos[3*(Seg[i].size-1)+1];
  z_pos[k] = Seg[i].Pos[3*(Seg[i].size-1)+2];
  r_sph[k] = Seg[i].Rvir[1];

  j = Seg[i].size+2;
  arvo_(&j,x_pos,y_pos,z_pos,r_sph,&vvv);

  Seg[i].vol  = Seg[i].Rvir[0]*Seg[i].Rvir[0]*Seg[i].Rvir[0];
  Seg[i].vol += Seg[i].Rvir[1]*Seg[i].Rvir[1]*Seg[i].Rvir[1];
  Seg[i].vol *= (4.0/3.0);
  Seg[i].vol *= M_PI;
  Seg[i].vol  = (float)vvv - Seg[i].vol;
  Seg[i].vol *= 1e-9;
  fac_len     = Seg[i].len - Seg[i].Rvir[0] - Seg[i].Rvir[1];
  fac_len    *= 1e-3;

  Seg[i].mass_part = cp.Mpart*(type_real)vec;
  Seg[i].mu        = Seg[i].mass_part/fac_len;
  Seg[i].rho       = Seg[i].mass_part/Seg[i].vol;

  free(x_pos);
  free(y_pos);
  free(z_pos);
  free(r_sph);

  return;

}
