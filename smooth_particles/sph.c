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

static type_int  **sum_fil;
static type_real **sum_mean_per;
static type_real **sum_quad_per;
static type_int *id_fil;
static type_int *k_fil;

static int compare_ids(const void *a, const void *b)
{
  //if((type_int *) a < (type_int *) b)
  //  return -1;

  //if((type_int *) a > (type_int *) b)
  //  return +1;

  //return 0;

  int *c, *d;
  
  c = (int *) a;
  d = (int *) b;
  
  return (*c - *d);
}

static void search_cent(const type_int idpar, const type_int tid, const type_real h2) 
{	
  type_int  nfil,id_cent_seg;
  type_int  *unique_fil;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox_grid[27], ibox;
  const type_real fac  = (type_real)grid.ngrid/cp.lbox;
  const type_real *Pos_cent = P[idpar].Pos;
  type_real r2, mm, Posprima[3], vprima[3], vdir[3];

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

  nfil = 0;
  ibox = 0;

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

        ibox_grid[ibox] = (ix * (long)grid.ngrid + iy) * (long)grid.ngrid + iz ;        
        nfil += grid.size[ibox_grid[ibox]];
        ibox++;
      }
    }  
  }

  unique_fil = (type_int *) malloc(nfil*sizeof(type_int));

  nfil = 0;
  for(ibox=0;ibox<27;ibox++)
  {
     id_cent_seg = grid.icell[ibox_grid[ibox]];

     while(id_cent_seg != grid.nobj)
     {
       unique_fil[nfil] = id_fil[id_cent_seg];
       id_cent_seg = grid.ll[id_cent_seg];
       nfil++;
     }
  }

  if(nfil>1)
  {
    qsort(unique_fil,nfil,sizeof(type_int),compare_ids);

    ix = 0;   
    for(iy=0; iy<nfil-1; iy++) 
      if(unique_fil[iy] != unique_fil[iy+1]) 
        unique_fil[ix++] = unique_fil[iy]; 
  
    unique_fil[ix++] = unique_fil[nfil-1];   
    nfil = ix;
  }

  for(ibox=0;ibox<nfil;ibox++)
  {
    Posprima[0] = Seg[unique_fil[ibox]].Pos_list[0] - Pos_cent[0];
    Posprima[1] = Seg[unique_fil[ibox]].Pos_list[1] - Pos_cent[1];
    Posprima[2] = Seg[unique_fil[ibox]].Pos_list[2] - Pos_cent[2];

    #ifdef PERIODIC
    if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
    if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
    if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
    if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
    if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
    if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
    #endif

    r2 = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

    if(r2 < Seg[unique_fil[ibox]].Rvir[0]) // uso el cuadrado   
      continue;

    Posprima[0] = Seg[unique_fil[ibox]].Pos_list[3*(Seg[unique_fil[ibox]].size)]   - Pos_cent[0];
    Posprima[1] = Seg[unique_fil[ibox]].Pos_list[3*(Seg[unique_fil[ibox]].size)+1] - Pos_cent[1];
    Posprima[2] = Seg[unique_fil[ibox]].Pos_list[3*(Seg[unique_fil[ibox]].size)+2] - Pos_cent[2];

    #ifdef PERIODIC
    if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
    if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
    if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
    if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
    if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
    if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
    #endif

    r2 = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

    if(r2 < Seg[unique_fil[ibox]].Rvir[1]) // uso el cuadrado   
      continue;

    mm = cp.lbox*cp.lbox;
    id_cent_seg = grid.nobj;

    for(ix=1;ix<Seg[unique_fil[ibox]].size;ix++)
    {

      r2 = 0.0f;
      vdir[0] = Seg[unique_fil[ibox]].Pos_list[3*ix]     - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)];
      vdir[1] = Seg[unique_fil[ibox]].Pos_list[3*ix + 1] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1) + 1];
      vdir[2] = Seg[unique_fil[ibox]].Pos_list[3*ix + 2] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1) + 2];
      #ifdef PERIODIC
      if(vdir[0]> 0.5*cp.lbox) vdir[0] -= cp.lbox;
      if(vdir[0]<-0.5*cp.lbox) vdir[0] += cp.lbox;
      if(vdir[1]> 0.5*cp.lbox) vdir[1] -= cp.lbox;
      if(vdir[1]<-0.5*cp.lbox) vdir[1] += cp.lbox;
      if(vdir[2]> 0.5*cp.lbox) vdir[2] -= cp.lbox;
      if(vdir[2]<-0.5*cp.lbox) vdir[2] += cp.lbox;
      #endif
      r2 = vdir[0]*vdir[0] + vdir[1]*vdir[1] + vdir[2]*vdir[2]; // tmb para mirar si caigo adentro del cilindro
      r2 = 1.0/sqrt(r2);

      vdir[0] *= r2; // directo
      vdir[1] *= r2; // directo
      vdir[2] *= r2; // directo

      Posprima[0] =  Pos_cent[0] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)]  ;
      Posprima[1] =  Pos_cent[1] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)+1];
      Posprima[2] =  Pos_cent[2] - Seg[unique_fil[ibox]].Pos_list[3*(ix-1)+2];

      #ifdef PERIODIC
      if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
      if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
      if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
      if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
      if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
      if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
      #endif

      vprima[1] = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2]; //r2

      if(vprima[1] <= h2 ) // entro en la esfera
      {
        vprima[0] = Posprima[0]*vdir[0] + Posprima[1]*vdir[1] + Posprima[2] * vdir[2]; // dot

        vprima[2] = fabs(vprima[1] - vprima[0]*vprima[0]);

        if(vprima[2]<mm)
        {
          mm = vprima[2];
          id_cent_seg = ix-1; // me quedo con el anterior
        }
      }
    }

    if(id_cent_seg == grid.nobj)
      continue;

    mm = r2 = 0.0;
    for(ix=0;ix<3;ix++)
    {
      vprima[ix] = P[idpar].Vel[ix]; // asigna

      //#ifdef CALCULA_VCM
        vprima[ix] -= Seg[unique_fil[ibox]].Vmedia[ix]; // asigna
      //#endif

      vdir[ix] = Seg[unique_fil[ibox]].Pos_list[3*(id_cent_seg+1)+ix] - Seg[unique_fil[ibox]].Pos_list[3*id_cent_seg+ix];
      #ifdef PERIODIC
      if(vdir[ix] >  0.5*cp.lbox) vdir[ix] -= cp.lbox;
      if(vdir[ix] < -0.5*cp.lbox) vdir[ix] += cp.lbox;
      #endif

      mm += (vprima[ix]*vprima[ix]);
      r2 += (vdir[ix]*vdir[ix]);
    }
    
    r2 = 1.0/sqrt(r2);

    vdir[0] *= r2;
    vdir[1] *= r2;
    vdir[2] *= r2;

    r2 = (vprima[0]*vdir[0] + vprima[1]*vdir[1] + vprima[2]*vdir[2]);
    mm -= (r2*r2);
    mm = sqrt(fabs(mm));

    sum_mean_per[tid][unique_fil[ibox]] += mm;
    sum_quad_per[tid][unique_fil[ibox]] += (mm*mm);  
    sum_fil[tid][unique_fil[ibox]]++;

    //mm = cp.lbox*cp.lbox;
    //id_cent_seg = grid.nobj;

    //for(ix=1;ix<(Seg[unique_fil[ibox]].size-1);ix++)
    //{
    //  Posprima[0] = Seg[unique_fil[ibox]].Pos_list[3*ix]   - Pos_cent[0];
    //  Posprima[1] = Seg[unique_fil[ibox]].Pos_list[3*ix+1] - Pos_cent[1];
    //  Posprima[2] = Seg[unique_fil[ibox]].Pos_list[3*ix+2] - Pos_cent[2];

    //  #ifdef PERIODIC
    //  if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
    //  if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
    //  if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
    //  if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
    //  if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
    //  if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
    //  #endif

    //  r2 = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

    //  if(r2<h2 && r2<mm)
    //  {
    //    mm = r2;
    //    id_cent_seg = ix;        
    //  }
    //}

    //if(id_cent_seg == grid.nobj)
    //  continue;
    //
    //// decide que vector usar

    //mm = r2 = 0.0;
    //for(ix=0;ix<3;ix++)
    //{
    //  Posprima[ix] = Pos_cent[ix] - Seg[unique_fil[ibox]].Pos_list[3*id_cent_seg+ix];
    //  #ifdef PERIODIC
    //  if(Posprima[ix] >  0.5*cp.lbox) Posprima[ix] -= cp.lbox;
    //  if(Posprima[ix] < -0.5*cp.lbox) Posprima[ix] += cp.lbox;
    //  #endif

    //  vdir[ix] = Seg[unique_fil[ibox]].Pos_list[3*(id_cent_seg+1)+ix] - Seg[unique_fil[ibox]].Pos_list[3*id_cent_seg+ix];
    //  #ifdef PERIODIC
    //  if(vdir[ix] >  0.5*cp.lbox) vdir[ix] -= cp.lbox;
    //  if(vdir[ix] < -0.5*cp.lbox) vdir[ix] += cp.lbox;
    //  #endif

    //  r2 += vdir[ix]*vdir[ix];
    //  mm += Posprima[ix]*vdir[ix]; // uso el signo del producto punto para saber a que versor se lo tengo que asignar
    //}

    //if(mm<0.0)
    //{
    //  id_cent_seg -= 1;

    //  r2 = 0.0;
    //  for(ix=0;ix<3;ix++)
    //  {

    //    Posprima[ix] = Pos_cent[ix] - Seg[unique_fil[ibox]].Pos_list[3*id_cent_seg+ix];
    //    #ifdef PERIODIC
    //    if(Posprima[ix] >  0.5*cp.lbox) Posprima[ix] -= cp.lbox;
    //    if(Posprima[ix] < -0.5*cp.lbox) Posprima[ix] += cp.lbox;
    //    #endif

    //    vdir[ix] = Seg[unique_fil[ibox]].Pos_list[3*(id_cent_seg+1)+ix] - Seg[unique_fil[ibox]].Pos_list[3*id_cent_seg+ix];
    //    #ifdef PERIODIC
    //    if(vdir[ix] >  0.5*cp.lbox) vdir[ix] -= cp.lbox;
    //    if(vdir[ix] < -0.5*cp.lbox) vdir[ix] += cp.lbox;
    //    #endif
    //    r2 += vdir[ix]*vdir[ix];
    //  }
    //}

    //r2 = 1.0/sqrt(r2);

    //vdir[0] *= r2;
    //vdir[1] *= r2;
    //vdir[2] *= r2;

    ////realizo el producto punto
    //mm = Posprima[0]*vdir[0] + Posprima[1]*vdir[1] + Posprima[2]*vdir[2];

    //r2 = 0.0f;
    //for(ix=0;ix<3;ix++)
    //{
    //  Posprima[ix] -= mm*vdir[ix];
    //
    //  #ifdef PERIODIC
    //  if(Posprima[ix]>  0.5f*cp.lbox) Posprima[ix] -= cp.lbox;
    //  if(Posprima[ix]< -0.5f*cp.lbox) Posprima[ix] += cp.lbox;
    //  #endif

    //  r2 += Posprima[ix]*Posprima[ix];
    //}
    //r2 = sqrt(r2);

    //mm = 0.0f;
    //for(ix=0;ix<3;ix++)
    //{
    //  vprima[ix] = P[idpar].Vel[ix]; // asigna

    //  //#ifdef CALCULA_VCM
    //    vprima[ix] -= Seg[unique_fil[ibox]].Vmedia[ix]; // asigna
    //  //#endif

    //  Posprima[ix] *= (1.0f/r2);
    //  mm += (vprima[ix]*Posprima[ix]);
    //}

    //sum_mean_per[tid][unique_fil[ibox]] += mm;
    //sum_quad_per[tid][unique_fil[ibox]] += (mm*mm);  
    //sum_fil[tid][unique_fil[ibox]]++;

  } //fin ibox

  free(unique_fil);
      
}

/*
static void search_cent(const type_real * __restrict__ Pos_cent, const type_int tid, const type_real h2) 
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

            sum_fil[tid][id_fil[id_cent_seg]]++;

          }// cierra el if
          
          id_cent_seg = grid.ll[id_cent_seg];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;

}
*/

extern void calculo_rho(void)
{
  long i, k, size;
  type_real aux_real, r[3];
  int k_arvo;

  aux_real  = cbrt(3.0/(4.0*M_PI*cp.Mpart*(type_real)cp.npart));
  aux_real *= cp.lbox;
  aux_real *= fof[1];
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

    Seg[i].Rvir[0] = RVIR_FACTOR*cbrt(Seg[i].Mass[0])*aux_real;
    Seg[i].Rvir[1] = RVIR_FACTOR*cbrt(Seg[i].Mass[1])*aux_real;
    //#ifdef CALCULA_VCM
      for(k=0;k<3;k++)
      {
        Seg[i].Vmedia[k] = (Seg[i].Mass[0]*Seg[i].Vnodos[k]+    \
                            Seg[i].Mass[1]*Seg[i].Vnodos[k+3])/ \
                           (Seg[i].Mass[0]+Seg[i].Mass[1]);
      }
    //#endif

    for(k=0;k<Seg[i].size;k++)
    {
      r[0] = Seg[i].Pos_list[3*k]   - Seg[i].Pos_list[0];
      r[1] = Seg[i].Pos_list[3*k+1] - Seg[i].Pos_list[1];
      r[2] = Seg[i].Pos_list[3*k+2] - Seg[i].Pos_list[2];

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
    k_arvo = Seg[i].size + 2;
    arvo_(&k_arvo, esferas, &vvv);
    Seg[i].vol = (float)vvv;

    //k = 0;
    esferas[0] = Seg[i].Pos_list[0];
    esferas[1] = Seg[i].Pos_list[1];
    esferas[2] = Seg[i].Pos_list[2];
    esferas[3] = Seg[i].Rvir[0];

    //k = 1;
    esferas[4] = Seg[i].Pos_list[3*(Seg[i].size-1)];
    esferas[5] = Seg[i].Pos_list[3*(Seg[i].size-1)+1];
    esferas[6] = Seg[i].Pos_list[3*(Seg[i].size-1)+2];
    esferas[7] = Seg[i].Rvir[1];

    k_arvo = 2; 
    vvv    = 0.0;
    arvo_(&k_arvo, esferas, &vvv);
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
  sum_mean_per = (type_real **) calloc(NTHREADS,sizeof(type_real *));
  sum_quad_per = (type_real **) calloc(NTHREADS,sizeof(type_real *));
  for(i=0;i<NTHREADS;i++)
  {
    sum_mean_per[i] = (type_real *) calloc(cp.nseg,sizeof(type_real));
    sum_quad_per[i] = (type_real *) calloc(cp.nseg,sizeof(type_real));
    sum_fil[i] = (type_int *) calloc(cp.nseg,sizeof(type_int));
  }

  grid.nobj = size;
  grid.ngrid = (long)(cp.lbox/(sqrt(2.0)*R_SPH));

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
  default(none) private(i,k) shared(cp,sum_fil,sum_mean_per,sum_quad_per,P,Seg,stdout) 
  for(i=0;i<cp.npart;i++)
  {
    if(i%10000000==0)
    {
      fprintf(stdout,"Particles %u %.4f\n",i,(float)i/(float)cp.npart);
      fflush(stdout);
    }

    k = omp_get_thread_num();
    search_cent(i,k,R_SPH*R_SPH);
  }

  #pragma omp parallel for num_threads(NTHREADS) default(none) \
  private(i,k,aux_real,r) shared(cp,Seg,sum_fil,sum_mean_per,sum_quad_per,stdout)
  for(i=0;i<cp.nseg;i++)
  {
    Seg[i].Rvir[0] = sqrt(Seg[i].Rvir[0]);
    Seg[i].Rvir[1] = sqrt(Seg[i].Rvir[1]);
    
    for(k=1;k<NTHREADS;k++)
    {
      sum_fil[0][i]         += sum_fil[k][i];
      sum_mean_per[0][i] += sum_mean_per[k][i];
      sum_quad_per[0][i] += sum_quad_per[k][i];
    }

    Seg[i].mass_part = cp.Mpart*(type_real)sum_fil[0][i];

    sum_mean_per[0][i] *= (sum_fil[0][i]>0) ? 1./(type_real)sum_fil[0][i] : 0.0f;

    Seg[i].sigma_per = sum_fil[0][i]>0 ? sqrt((sum_quad_per[0][i] - \
    sum_fil[0][i]*sum_mean_per[0][i]*sum_mean_per[0][i])*(1.f/(type_real)(sum_fil[0][i]-1))) : 0.0f;

    aux_real = Seg[i].len - Seg[i].Rvir[0] - Seg[i].Rvir[1];
    aux_real *= 1e-3;

    Seg[i].mu  = aux_real<0.0 ? 0.0f : Seg[i].mass_part/aux_real;
    Seg[i].rho = aux_real<0.0 ? 0.0f : Seg[i].mass_part/Seg[i].vol;

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

  for(i=0;i<NTHREADS;i++)
  {
    free(sum_fil[i]);
    free(sum_mean_per[i]);
    free(sum_quad_per[i]);
  }
  free(sum_fil);
  free(sum_mean_per);
  free(sum_quad_per);

  free(id_fil);
  free(k_fil);
  grid_free();
  free(P);

  return;
}
