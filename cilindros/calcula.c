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
#include "calcula.h"
#include "colores.h"
#include "leesnap.h"
#ifdef SAVEPART
  #include "bitmask.h"
#endif

static void set_name(char * name, const type_int NNN, char * prefix)
{

  #ifdef FIXED_SEPARATION
    sprintf(name,"%.4u_fixed",NNN);
  #else
    sprintf(name,"%.4u_varia",NNN);
  #endif

  #ifdef MCRITIC
    sprintf(name,"%s_%s_cut",name,prefix);
  #else
    sprintf(name,"%s_%s",name,prefix);
  #endif

  #ifdef BIN_LOG
    sprintf(name,"%s_LOG",name);
  #else
    sprintf(name,"%s_LIN",name);
  #endif

  return;

}

static inline type_int point_inside(type_real dot, type_real rlong_2)
{
  if(fabs(dot)>rlong_2){

      return 0;

  }else{

    return 1;

  }
}

static int cmp(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

static void stadistic(type_int n, type_real *MAX, type_real *MIN, type_real *LMAX) 
{
  type_int i, j;
  type_real mean, median, sumquad, std, *a;

  a = (type_real *) malloc(n*sizeof(type_real));

  j = 0;
  mean = sumquad = 0.0;
  for(i=0;i<n;i++) 
  {
    a[j] = Seg[i].len;
    mean+=a[j];
    sumquad+=a[j]*a[j];
    j++;
  }

  if(j!=n)
     a = (type_real *) realloc(a,j*sizeof(type_real));

  qsort(a, j, sizeof(type_real), cmp);

  mean /= (type_real)j;
  sumquad /= (type_real)j;
  
  std = sqrt(sumquad - mean*mean);
  
  fprintf(stdout,"elements %d\n",j);
  fprintf(stdout,"mean %f\n",mean);
  fprintf(stdout,"std  %f\n",std);
  fprintf(stdout,"min %f\n",a[0]);
  fprintf(stdout,"max %f\n",a[j-1]);

  *LMAX = a[j-1]/(type_real)nbins;
  *LMAX *= 2.0;

  if(j%2==0) {
      median = (a[j/2] + a[j/2 - 1])*0.5;
      fprintf(stdout,"median %f\n",median);
  } else {
      median = a[j/2];
      fprintf(stdout,"median %f\n",median);
  }

  if(mean>median)
    *MAX = mean;
  else
    *MAX = median;

  //*MAX = 500.*ceil(*MAX/1000.);
  *MAX = 40000.;  
  *MIN = *MAX/50.;

  free(a);
  
  return;
}

#ifdef BIN_LOG

  static void logspace(type_real *rcil2, const type_real max, const type_real min, const type_int bins) 
  {
  
    type_int i;
    type_real end   = log10(max);
    type_real start = log10(min);
    type_real delta = (end - start) / (bins-1);
  
    rcil2[0] = pow(10.0,end);   
    rcil2[0] *= rcil2[0];
  
    for(i=1;i<bins; i++)
    {
      rcil2[i] = pow(10.0,(start + delta * (bins - 1 - i)));
      rcil2[i] *= rcil2[i];
    }  
  
    return;
  }

#else

  static void linspace(type_real *rcil2, const type_real max, const type_real min, const type_int bins) 
  {
  
    type_int i;
    type_real end   = max;
    type_real start = min;
    type_real delta = (end - start) / (bins-1);
  
    rcil2[0] = end;   
    rcil2[0] *= rcil2[0];
  
    for(i=1;i<bins; i++)
    {
      rcil2[i] = start + delta * (bins - 1 - i);
      rcil2[i] *= rcil2[i];
    }  
  
    return;
  }

#endif


#ifdef SAVEPART
  static static void calc_part(type_real * const Pos_cent, type_real * const Vmedia, type_int * vec, type_int * ncont, type_int * npart, type_real ** mean, type_real ** quad, \
  type_real * versor, type_real * mean_par, type_real * quad_par, type_real * mean_perp, type_real * quad_perp, type_real * const rcil2, type_real rlong_2, type_int binsper)
#else        
  static void calc_part(type_real * const Pos_cent, type_real * const Vmedia, type_int * npart, type_real ** mean, type_real ** quad, \
  type_real * versor, type_real * mean_par, type_real * quad_par, type_real * mean_perp, type_real * quad_perp, type_real * const rcil2, type_real rlong_2, type_int binsper)
#endif
{
  type_int i, j, k;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox, ngrid, ilim;
  type_real lbox,fac,lbox2;
  type_real dot,dis;
  type_real aux_real, aux_dot, aux_mod;
  type_real xxx[3];
  type_real Posprima[3],Vprima[3];
  #ifdef HUECOS
  type_int bin;
  #endif

  ////////////////////////////////////////////////////
  memset(npart,0,binsper*sizeof(type_int));
  memset(mean_par,0.0,binsper*sizeof(type_real));
  memset(quad_par,0.0,binsper*sizeof(type_real));
  memset(mean_perp,0.0,binsper*sizeof(type_real));
  memset(quad_perp,0.0,binsper*sizeof(type_real));
  for(j=0;j<binsper;j++)
  {
    memset(mean[j],0.0,3*sizeof(type_real));
    memset(quad[j],0.0,3*sizeof(type_real));
  }
  ////////////////////////////////////////////////////

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  dis = sqrt(rcil2[0]);
  ilim = dis>rlong_2 ? (long)(dis*fac)+1 : (long)(rlong_2*fac)+1;

  ixc  = (long)(Pos_cent[0]*fac);
  ixci = ixc - ilim;
  ixcf = ixc + ilim;
  iyc  = (long)(Pos_cent[1]*fac);
  iyci = iyc - ilim;
  iycf = iyc + ilim;
  izc  = (long)(Pos_cent[2]*fac);
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

        #ifndef REORDER
          i = grid.llirst[ibox];
        #endif
      
        #ifdef REORDER
          for(i=grid.icell[ibox];i<grid.icell[ibox]+grid.size[ibox];i++)
        #else
          while(i != cp.npart)
        #endif
        {
          #ifdef COLUMN
            Posprima[0] = P.x[i] - Pos_cent[0];
            Posprima[1] = P.y[i] - Pos_cent[1];
            Posprima[2] = P.z[i] - Pos_cent[2];
          #else
            Posprima[0] = P[i].Pos[0] - Pos_cent[0];
            Posprima[1] = P[i].Pos[1] - Pos_cent[1];
            Posprima[2] = P[i].Pos[2] - Pos_cent[2];
          #endif

          #ifdef PERIODIC
          if(Posprima[0] >  lbox2) Posprima[0] = Posprima[0] - lbox;
          if(Posprima[1] >  lbox2) Posprima[1] = Posprima[1] - lbox;
          if(Posprima[2] >  lbox2) Posprima[2] = Posprima[2] - lbox;
          if(Posprima[0] < -lbox2) Posprima[0] = Posprima[0] + lbox;
          if(Posprima[1] < -lbox2) Posprima[1] = Posprima[1] + lbox;
          if(Posprima[2] < -lbox2) Posprima[2] = Posprima[2] + lbox;
          #endif

          dot = Posprima[0]*versor[0]+Posprima[1]*versor[1]+Posprima[2]*versor[2];

          if(point_inside(dot,rlong_2)==1)
          {

            dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2] - dot*dot;

            #ifdef HUECOS
            //if(dis<rcil2[0] && dis>rcil2[binsper])
            if(dis<rcil2[0])
            {

              bin = 0;
              for(j=1;j<binsper;j++)
              {
                bin = dis<rcil2[j] ? j : bin;
              }

              #ifdef SAVEPART
              if(bin>=3 && !TestBit(vec,i))
              {
                SetBit(vec,i);
                *ncont = *ncont + 1;
              }
              #endif

              #ifdef COLUMN
              Vprima[0] = P.vx[i] - Vmedia[3*bin+0];
              Vprima[1] = P.vy[i] - Vmedia[3*bin+1];
              Vprima[2] = P.vz[i] - Vmedia[3*bin+2];
              #else
              Vprima[0] = P[i].Vel[0] - Vmedia[3*bin+0];
              Vprima[1] = P[i].Vel[1] - Vmedia[3*bin+1];
              Vprima[2] = P[i].Vel[2] - Vmedia[3*bin+2];
              #endif

              aux_mod = aux_dot = 0.0;
              for(k=0;k<3;k++)
              {
                mean[bin][k] += Vprima[k];
                quad[bin][k] += (Vprima[k]*Vprima[k]);
  
                aux_dot += (Vprima[k]*versor[k]);
  
                xxx[k] = (Posprima[k] - dot*versor[k]);                               
  
                #ifdef PERIODIC
                if(xxx[k]>  lbox2) xxx[k] -= lbox;
                if(xxx[k]< -lbox2) xxx[k] += lbox;
                #endif
  
                aux_mod += (xxx[k]*xxx[k]);
              }
  
              mean_par[bin] += aux_dot;
              quad_par[bin] += (aux_dot*aux_dot);
              aux_mod = 1.0/sqrt(aux_mod);
  
              aux_real = 0.0;
              for(k=0;k<3;k++)
              {
                xxx[k] *= aux_mod;
                aux_real += (Vprima[k]*xxx[k]);
              }              
  
              mean_perp[bin] += aux_real;
              quad_perp[bin] += (aux_real*aux_real);
   
              npart[bin]++;
  
            }

            #else

            for(j=0;j<binsper;j++)
            {
              if(dis<rcil2[j])
              {
                if(j==0)
                {

                  #ifdef SAVEPART

                  if(!TestBit(vec,i))
                  {
                    SetBit(vec,i);
                    *ncont = *ncont + 1;
                  }
                  #endif

                  aux_mod = aux_dot = 0.0;
                  for(k=0;k<3;k++)
                  {
                    Vprima[k] = P[i].Vel[k] - Vmedia[j+binsper*k];

                    mean[j][k] += Vprima[k];
                    quad[j][k] += (Vprima[k]*Vprima[k]);

                    aux_dot += (Vprima[k]*versor[k]);

                    xxx[k] = (Posprima[k] - dot*versor[k]);                               

                    #ifdef PERIODIC
                    if(xxx[k]>  lbox2) xxx[k] -= lbox;
                    if(xxx[k]< -lbox2) xxx[k] += lbox;
                    #endif

                    aux_mod += (xxx[k]*xxx[k]);
                  }

                  mean_par[j] += aux_dot;
                  quad_par[j] += (aux_dot*aux_dot);
                  aux_mod = 1.0/sqrt(aux_mod);

                  aux_real = 0.0;
                  for(k=0;k<3;k++)
                  {
                    xxx[k] *= aux_mod;
                    aux_real += (Vprima[k]*xxx[k]);
                  }              

                  mean_perp[j] += aux_real;
                  quad_perp[j] += (aux_real*aux_real);
 
                }else{

                  for(k=0;k<3;k++)
                  {
                    mean[j][k] += Vprima[k];
                    quad[j][k] += (Vprima[k]*Vprima[k]);
                  }              

                  mean_par[j] += aux_dot;
                  quad_par[j] += (aux_dot*aux_dot);

                  mean_perp[j] += aux_real;
                  quad_perp[j] += (aux_real*aux_real); 
                  
                }

                npart[j]++;

              }else{

                break;

              }
            }

          #endif

          }// cierra el if

          #ifndef REORDER
            i = grid.ll[i];
          #endif

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

#ifdef SAVEPART

  static void calcular_mean(const type_int i, const type_int binsper, const type_real rsep, \
           const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
  type_real **vmean, type_real **vquad, type_real *vmean_par, type_real *vquad_par,\
  type_real *vmean_perp, type_real *vquad_perp, type_real *nodo, struct node_sph **root,\
  type_int *size_part, type_int *vpart)

#else

  static void calcular_mean(const type_int i, const type_int binsper, const type_real rsep, \
           const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
  type_real **vmean, type_real **vquad, type_real *vmean_par, type_real *vquad_par,\
  type_real *vmean_perp, type_real *vquad_perp, type_real *nodo, struct node_sph **root)

#endif
{
  type_int k,l,j,idim;
  type_int bin,lbin;
  type_int vsize[ncil];
  type_real vdir[3];
  type_real Pos_cent[3];
  type_real mod_v[3];
  type_real r,rbin,racum;
  type_real *volinv;                // En Kpc
  type_real dens, lenrbin;

  volinv = (type_real *) malloc(ncil*sizeof(type_real));
  for(j=0;j<ncil;j++)
  {
    volinv[j] = M_PI*rlong*rcil2[j];
    
    #ifdef HUECOS
      r = M_PI*rlong*rcil2[j+1];
      if(j != ncil-1) // para hacer el ultimo cilindro entero
        volinv[j] -= r;
    #endif  

    volinv[j] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[j])/(type_real)cp.npart;
  }

  racum   = 0.0;
  lenrbin = 0.0;
  l = 0;

  ///////////////////////////////////////////////////////
  r = 0.0;
  for(idim=0;idim<3;idim++)
  {
    vdir[idim] = Gr[Seg[i].list[1]].Pos[idim]-Gr[Seg[i].list[0]].Pos[idim];    
    #ifdef PERIODIC
    if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
    if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
    #endif
    r += vdir[idim]*vdir[idim];
  }
  r = sqrt(r);

  for(idim=0;idim<3;idim++)
    vdir[idim] *= (1./r);
  ///////////////////////////////////////////////////////

  #ifdef SAVEPART
  calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
  #else
  calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
  #endif

  nodo[l] = lenrbin; 

  for(j=0;j<ncil;j++)
  {
    dens = (type_real)vsize[j]*volinv[j];
    push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
    vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
  }
  l++; // sumo

  for(k=1; k<Seg[i].size; k++)
  {    
    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[i].list[k]].Pos[idim]-Gr[Seg[i].list[k-1]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);

    for(idim=0;idim<3;idim++)
    {
      vdir[idim] *= (1.0/r);
      Pos_cent[idim] = Gr[Seg[i].list[k-1]].Pos[idim];
    }

    rbin = rsep-racum;

    if(r<rbin)
    {
      racum += r;
      continue;
    }
  
    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] += rbin*vdir[idim];
      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

    #ifdef SAVEPART
    calc_part(Pos_cent,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #else
    calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #endif
    lenrbin += rsep;

    nodo[l] = lenrbin; 
    for(j=0;j<ncil;j++)
    {
      dens = (type_real)vsize[j]*volinv[j];
      push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
      vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
    }
    l++; // sumo

    r = 0.0;
    for(idim=0;idim<3;idim++)
    {
      mod_v[idim] = Gr[Seg[i].list[k]].Pos[idim]-Pos_cent[idim];
      #ifdef PERIODIC
      mod_v[idim] = mod_v[idim] >= 0.5*cp.lbox ? mod_v[idim]-cp.lbox : mod_v[idim];
      mod_v[idim] = mod_v[idim] < -0.5*cp.lbox ? mod_v[idim]+cp.lbox : mod_v[idim];
      #endif
      r += mod_v[idim]*mod_v[idim];
    }
    r = sqrt(r);

    lbin = (type_int)(r/rsep);

    for(bin=0;bin<lbin;bin++)
    {

      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }

      #ifdef SAVEPART
      calc_part(Pos_cent,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #else
      calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #endif

      lenrbin += rsep;
      nodo[l] = lenrbin; 

      for(j=0;j<ncil;j++)
      {
        dens = (type_real)vsize[j]*volinv[j];
        push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
        vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
      }

      l++; // sumo
    }

    racum = r-(type_real)lbin*rsep;
  }

  if(l==nbins) 
    l--; // sobre escribo el ultimo

  ///////////////////////////////////////////////////////
  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    vdir[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim]-Gr[Seg[i].list[Seg[i].size-2]].Pos[idim];    
    #ifdef PERIODIC
    if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
    if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
    #endif
    r += vdir[idim]*vdir[idim];
  }
  r = sqrt(r);

  for(idim=0;idim<3;idim++)
    vdir[idim] *= (1./r);
  ///////////////////////////////////////////////////////

  #ifdef SAVEPART
  calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
  #else
  calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
  #endif

  nodo[l] = Seg[i].len; 
  for(j=0;j<ncil;j++)
  {
    dens = (type_real)vsize[j]*volinv[j];
    push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
    vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
  }
  l++; // sumo

  assert(l==nbins);

  free(volinv);
}

#ifdef CALCULA_MEDIA

static void calc_media(type_real * const Pos_cent, type_real * const versor, type_int *numpart, type_real *vel_media, \
                       type_real * const rcil2, const type_real rlong_2, const type_int binsper)
{
  type_int i,j,bin;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox, ngrid, ilim;
  type_real Posprima[3];
  type_real lbox,fac,lbox2;
  type_real dot,dis;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  dis = sqrt(rcil2[0]);
  ilim = dis>rlong_2 ? (long)(dis*fac)+1 : (long)(rlong_2*fac)+1;

  ixc  = (long)(Pos_cent[0]*fac);
  ixci = ixc - ilim;
  ixcf = ixc + ilim;
  iyc  = (long)(Pos_cent[1]*fac);
  iyci = iyc - ilim;
  iycf = iyc + ilim;
  izc  = (long)(Pos_cent[2]*fac);
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

        #ifndef REORDER
          i = grid.llirst[ibox];
        #endif
      
        #ifdef REORDER
          for(i=grid.icell[ibox];i<grid.icell[ibox]+grid.size[ibox];i++)
        #else
          while(i != cp.npart)
        #endif
        {
          #ifdef COLUMN
            Posprima[0] = P.x[i] - Pos_cent[0];
            Posprima[1] = P.y[i] - Pos_cent[1];
            Posprima[2] = P.z[i] - Pos_cent[2];
          #else
            Posprima[0] = P[i].Pos[0] - Pos_cent[0];
            Posprima[1] = P[i].Pos[1] - Pos_cent[1];
            Posprima[2] = P[i].Pos[2] - Pos_cent[2];
          #endif

          #ifdef PERIODIC
          if(Posprima[0] >  lbox2) Posprima[0] = Posprima[0] - lbox;
          if(Posprima[1] >  lbox2) Posprima[1] = Posprima[1] - lbox;
          if(Posprima[2] >  lbox2) Posprima[2] = Posprima[2] - lbox;
          if(Posprima[0] < -lbox2) Posprima[0] = Posprima[0] + lbox;
          if(Posprima[1] < -lbox2) Posprima[1] = Posprima[1] + lbox;
          if(Posprima[2] < -lbox2) Posprima[2] = Posprima[2] + lbox;
          #endif

          dot = Posprima[0]*versor[0]+Posprima[1]*versor[1]+Posprima[2]*versor[2];

          if(point_inside(dot,rlong_2)==1)
          {
            dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2] - dot*dot;

            //if(dis<rcil2[0] && dis>rcil2[binsper])
            if(dis<rcil2[0])
            {
              bin = 0;
              for(j=1;j<binsper;j++)
              {
                bin = dis<rcil2[j] ? j : bin;
              }

              #ifdef COLUMN
              vel_media[3*bin+0] += P.vx[i];
              vel_media[3*bin+1] += P.vy[i];
              vel_media[3*bin+2] += P.vz[i];
              #else
              vel_media[3*bin+0] += P[i].Vel[0];
              vel_media[3*bin+1] += P[i].Vel[1];
              vel_media[3*bin+2] += P[i].Vel[2];
              #endif
              numpart[bin]++;

            }
          } // cierra dis

          #ifndef REORDER
            i = grid.ll[i];
          #endif

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

static void cylinder_mean(const type_int i, const type_real rsep, const type_real rlong, \
                   const type_real rlong_2, type_real * const rcil2)
{
  type_int *npart;
  type_int k,l,idim;
  type_int bin,lbin;
  type_real vdir[3];
  type_real Pos_cent[3];
  type_real mod_v[3];
  type_real r,rbin,racum;

  npart = (type_int *) calloc(ncil,sizeof(type_int));
  l = 0;
  racum = 0.0f;

  ///////////////////////////////////////////////////////
  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    vdir[idim] = Gr[Seg[i].list[1]].Pos[idim]-Gr[Seg[i].list[0]].Pos[idim];    
    #ifdef PERIODIC
    if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
    if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
    #endif
    r += vdir[idim]*vdir[idim];
  }
  r = sqrt(r);

  for(idim=0;idim<3;idim++)
    vdir[idim] *= (1.0/r);
  ///////////////////////////////////////////////////////

  calc_media(Gr[Seg[i].list[0]].Pos,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);             
  l++;

  for(k=1; k<Seg[i].size; k++)
  {

    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[i].list[k]].Pos[idim]-Gr[Seg[i].list[k-1]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);

    for(idim=0;idim<3;idim++)
    {
      vdir[idim] *= (1.0/r);
      Pos_cent[idim] = Gr[Seg[i].list[k-1]].Pos[idim];
    }

    rbin = rsep-racum;

    if(r<rbin)
    {
      racum += r;
      continue;
    }

    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] += rbin*vdir[idim];  
      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

    calc_media(Pos_cent,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);             
    l++;

    r = 0.0;
    for(idim=0;idim<3;idim++)
    {
      mod_v[idim] = Gr[Seg[i].list[k]].Pos[idim]-Pos_cent[idim];
      #ifdef PERIODIC
      mod_v[idim] = mod_v[idim] >= 0.5*cp.lbox ? mod_v[idim]-cp.lbox : mod_v[idim];
      mod_v[idim] = mod_v[idim] < -0.5*cp.lbox ? mod_v[idim]+cp.lbox : mod_v[idim];
      #endif
      r += mod_v[idim]*mod_v[idim];
    }
    r = sqrt(r);

    lbin = (type_int)(r/rsep);

    for(bin=0;bin<lbin;bin++)
    {

      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];  
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }
      calc_media(Pos_cent,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);

      l++;
    }

    racum = r-(type_real)lbin*rsep;
  }

  if(l==nbins)
    l--;     // sobre escribo el ultimo

  ///////////////////////////////////////////////////////
  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    vdir[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim]-Gr[Seg[i].list[Seg[i].size-2]].Pos[idim];    
    #ifdef PERIODIC
    if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
    if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
    #endif
    r += vdir[idim]*vdir[idim];
  }
  r = sqrt(r);

  for(idim=0;idim<3;idim++)
    vdir[idim] *= (1./r);
  ///////////////////////////////////////////////////////

  calc_media(Gr[Seg[i].list[Seg[i].size-1]].Pos,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);             
  l++;

  assert(l==nbins);

  for(bin=0;bin<ncil;bin++)
  {
    Seg[i].Vmedia[3*bin + 0] *= (npart[bin]>0 ? 1.0f/(type_real)npart[bin] : 0.0);
    Seg[i].Vmedia[3*bin + 1] *= (npart[bin]>0 ? 1.0f/(type_real)npart[bin] : 0.0);
    Seg[i].Vmedia[3*bin + 2] *= (npart[bin]>0 ? 1.0f/(type_real)npart[bin] : 0.0);
  }

  free(npart);

}

#endif

#ifdef EXTEND

  #ifdef SAVEPART
  
    static void ext_mean(const type_int i, const type_int binsper, const type_real rsep,\
    const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
    type_int *size_part, type_int *vpart, struct dat_struct * dat)

  #else
  
    static void ext_mean(const type_int i, const type_int binsper, const type_real rsep,\
    const type_real rlong, const type_real rlong_2, type_real * const rcil2, struct dat_struct * dat)

  #endif
  {
    type_int l,j,idim,lbin;
    type_int vsize[ncil];
    type_int * npart;
    type_real aux_mean, aux_rms;
    type_real r, dens;
    type_real vdir[3];
    type_real Pos_cent[3];
    type_real *volinv;                // En Kpc
    type_real *nodo;
    type_real *vmean_par, *vquad_par;
    type_real *vmean_perp, *vquad_perp;
    type_real **vmean, **vquad;
    struct node_sph **root;
  
    ////////////////////////////////////////////////////////////////////////////////
    l = 0;
    lbin = nbins%2==0 ? nbins/2 : (type_int)floor((float)nbins/2);
    assert(lbin!=0);

    volinv = (type_real *) malloc(ncil*sizeof(type_real));
    npart  = (type_int *) calloc(ncil,sizeof(type_int));

    nodo = (type_real *) malloc(lbin*sizeof(type_real));
    vmean_par = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_par = (type_real *) malloc(ncil*sizeof(type_real));
    vmean_perp = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_perp = (type_real *) malloc(ncil*sizeof(type_real));   

    root = (struct node_sph **) malloc(ncil*sizeof(struct node_sph *));
    vmean = (type_real **) malloc(ncil*sizeof(type_real *));
    vquad = (type_real **) malloc(ncil*sizeof(type_real *));
    for(j=0;j<ncil;j++)
    {
      root[j] = (struct node_sph *) malloc(lbin*sizeof(struct node_sph));
      vmean[j] = (type_real *) malloc(3*sizeof(type_real));
      vquad[j] = (type_real *) malloc(3*sizeof(type_real));
    }
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    for(j=0;j<ncil;j++) 
    {
      volinv[j] = M_PI*rlong*rcil2[j];
  
      #ifdef HUECOS
        r = M_PI*rlong*rcil2[j+1];
        if(j != ncil-1) // para hacer el ultimo cilindro entero
          volinv[j] -= r;
      #endif  
  
      volinv[j] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[j])/(type_real)cp.npart;
    }
    ////////////////////////////////////////////////////////////////////////////////

    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[i].list[1]].Pos[idim]-Gr[Seg[i].list[0]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);

    l = 0;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] *= (1/r); // directo
      Pos_cent[idim] = Gr[Seg[i].list[0]].Pos[idim] - (lbin-1)*rsep*vdir[idim]; // cambio el origen
      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

    #ifdef SAVEPART
      calc_part(Pos_cent,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
    #else
      calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #endif

    nodo[l] = -(type_real)(lbin-1)*rsep;
    for(j=0;j<ncil;j++)
    {
      dens = (type_real)vsize[j]*volinv[j];
      push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
      vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
    }
    l++; // sumo  
 
    while(l<lbin)
    {  
      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }

      #ifdef SAVEPART
        if(l==lbin-1)
          calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
        else
          calc_part(Pos_cent,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #else
        if(l==lbin-1)
          calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
        else
          calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #endif
 
      nodo[l] = nodo[l-1]+rsep;   
      for(j=0;j<ncil;j++)
      {
        dens = (type_real)vsize[j]*volinv[j];
        push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
        vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
      }
  
      l++; // sumo
    }
  
    assert(l==lbin);
    
    ////////////////////////////////////////////////////////////////////////////////
    for(j=0;j<lbin;j++)
    {
      nodo[j] = j==lbin-1 ? 0.0f : nodo[j];
      //r = nodo[j]/Seg[i].len; // Para normalizarlo
      r = nodo[j];              // Sin normalizacion
      dat[i].nodo[nbins+j] = nodo[j];
    }

    for(j=0;j<ncil;j++)
    {
      r  = sqrt(rcil2[j]);
      r += sqrt(rcil2[j+1]);
      r *= 0.5;

      for(l=0;l<lbin;l++)
      {
        root[j][l].b -= 1.0;

        dat[i].npart[j][nbins+l]     = root[j][l].a; // npart
        dat[i].rho_delta[j][nbins+l] = root[j][l].b; // rho - delta

        ////////////////////////////////////////////////////////////////////////////////////

        for(idim=0;idim<3;idim++)
        {
          root[j][l].c[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        }
        root[j][l].e *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].g *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;

        aux_mean = aux_rms = 0.0;

        for(idim=0;idim<3;idim++)
        {
          aux_mean  += root[j][l].c[idim];
          aux_rms   += root[j][l].d[idim] - (type_real)root[j][l].a*root[j][l].c[idim]*root[j][l].c[idim];
        }
        
        aux_mean /= 3.0;        
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms)*(1.f/(type_real)(root[j][l].a-1)))/3.0 : 0.0;
        //aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms))/3.0 : 0.0;

        dat[i].mean[j][nbins+l] = aux_mean; // mean
        dat[i].rms[j][nbins+l]  = aux_rms;  // rms

        aux_mean = root[j][l].e;
        aux_rms  = root[j][l].a>10 ? sqrt((fabs(root[j][l].f - \
        (type_real)root[j][l].a*root[j][l].e*root[j][l].e))*(1.f/(type_real)(root[j][l].a-1))) : 0.0;

        dat[i].mean_par[j][nbins+l] = aux_mean; // mean par 
        dat[i].rms_par[j][nbins+l]  = aux_rms;  // rms  par

        aux_mean = root[j][l].g;
        aux_rms  = root[j][l].a>10 ? sqrt((fabs(root[j][l].h - \
        (type_real)root[j][l].a*root[j][l].g*root[j][l].g))*(1.f/(type_real)(root[j][l].a-1))) : 0.0;

        dat[i].mean_perp[j][nbins+l] = aux_mean;// mean perp        
        dat[i].rms_perp[j][nbins+l]  = aux_rms; // rms  perp
      }

      free(vmean[j]);
      free(vquad[j]);
      free(root[j]);
    }

    free(nodo);
    free(root);
    free(vmean);
    free(vquad);
    free(vmean_par);
    free(vquad_par);
    free(vmean_perp);
    free(vquad_perp);
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    l = 0;
    //memset(Vmean,0.0,3*ncil*sizeof(type_real));
    memset(npart,0,ncil*sizeof(type_int));
    nodo = (type_real *) malloc(lbin*sizeof(type_real));
    vmean_par = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_par = (type_real *) malloc(ncil*sizeof(type_real));
    vmean_perp = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_perp = (type_real *) malloc(ncil*sizeof(type_real));   

    root = (struct node_sph **) malloc(ncil*sizeof(struct node_sph *));
    vmean = (type_real **) malloc(ncil*sizeof(type_real *));
    vquad = (type_real **) malloc(ncil*sizeof(type_real *));
    for(j=0;j<ncil;j++)
    {
      root[j] = (struct node_sph *) malloc(lbin*sizeof(struct node_sph));
      vmean[j] = (type_real *) malloc(3*sizeof(type_real));
      vquad[j] = (type_real *) malloc(3*sizeof(type_real));
    }
    ////////////////////////////////////////////////////////////////////////////////

    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim]-Gr[Seg[i].list[Seg[i].size-2]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);

    for(idim=0;idim<3;idim++)
      vdir[idim] *= (1/r); // directo

    l = 0;
    for(idim=0;idim<3;idim++)
      Pos_cent[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim];

    #ifdef SAVEPART
    calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
    #else
    calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #endif
  
    nodo[l] = 0.0f;   
    for(j=0;j<ncil;j++)
    {
      dens = (type_real)vsize[j]*volinv[j];
      push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
      vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
    }
    l++; // sumo  
 
    while(l<lbin)
    {  
      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }
  
      #ifdef SAVEPART
        calc_part(Pos_cent,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #else
        calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #endif
  
      nodo[l] = nodo[l-1]+rsep;   
      for(j=0;j<ncil;j++)
      {
        dens = (type_real)vsize[j]*volinv[j];
        push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
        vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
      }
  
      l++; // sumo
    }
  
    assert(l==lbin);
   
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(j=0;j<lbin;j++)
    {
      nodo[j] = j==0 ? Seg[i].len : Seg[i].len+nodo[j];
      //r = nodo[j]/Seg[i].len; // Para normalizarlo
      r = nodo[j];              // Sin normalizacion
      //fwrite(&r,sizeof(type_real),1,pfextend);
      dat[i].nodo[nbins+lbin+j] = nodo[j];
    }

    for(j=0;j<ncil;j++)
    {
      r  = sqrt(rcil2[j]);
      r += sqrt(rcil2[j+1]);
      r *= 0.5;

      for(l=0;l<lbin;l++)
      {
        root[j][l].b -= 1.0;

        dat[i].npart[j][nbins+lbin+l]     = root[j][l].a; // npart
        dat[i].rho_delta[j][nbins+lbin+l] = root[j][l].b; // rho - delta

        ////////////////////////////////////////////////////////////////////////////////////

        for(idim=0;idim<3;idim++)
        {
          root[j][l].c[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        }
        root[j][l].e *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].g *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;

        aux_mean = aux_rms = 0.0;

        for(idim=0;idim<3;idim++)
        {
          aux_mean  += root[j][l].c[idim];
          aux_rms   += root[j][l].d[idim] - (type_real)root[j][l].a*root[j][l].c[idim]*root[j][l].c[idim];
        }
        
        aux_mean /= 3.0;        
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms)*(1.f/(type_real)(root[j][l].a-1)))/3.0 : 0.0;
        //aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms))/3.0 : 0.0;

        dat[i].mean[j][nbins+lbin+l] = aux_mean; // mean
        dat[i].rms[j][nbins+lbin+l]  = aux_rms;  // rms

        aux_mean = root[j][l].e;
        aux_rms  = root[j][l].a>10 ? sqrt((fabs(root[j][l].f - \
        (type_real)root[j][l].a*root[j][l].e*root[j][l].e))*(1.f/(type_real)(root[j][l].a-1))) : 0.0;

        dat[i].mean_par[j][nbins+lbin+l] = aux_mean; // mean par 
        dat[i].rms_par[j][nbins+lbin+l]  = aux_rms;  // rms  par

        aux_mean = root[j][l].g;
        aux_rms  = root[j][l].a>10 ? sqrt((fabs(root[j][l].h - \
        (type_real)root[j][l].a*root[j][l].g*root[j][l].g))*(1.f/(type_real)(root[j][l].a-1))) : 0.0;

        dat[i].mean_perp[j][nbins+lbin+l] = aux_mean;// mean perp        
        dat[i].rms_perp[j][nbins+lbin+l]  = aux_rms; // rms  perp
      }

      free(vmean[j]);
      free(vquad[j]);
      free(root[j]);
    }

    free(nodo);
    free(root);
    free(vmean);
    free(vquad);
    free(vmean_par);
    free(vquad_par);
    free(vmean_perp);
    free(vquad_perp);
    free(volinv);
    free(npart);
  }

#endif

#ifdef VEL_RELATIVA

static void relativa_write(const type_int NNN, type_real *fof)
{
  type_int i, idim;
  type_real vdir[3];
  type_real r;
  FILE *pfrela;
  char filename[200], name[200];

  BLUE("******************************\n");

  sprintf(filename,"Escribe el archivo para calcular las velocidades relativas\n");GREEN(filename);

  set_name(name,NNN,"vrelativa");
  sprintf(filename,"../%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pfrela = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(int),1,pfrela);        
 
  ///////////////////////////////////////////////////////
  
  for(i=0;i<cp.nseg;i++)
  {
    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim]-Gr[Seg[i].list[0]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);

    for(idim=0;idim<3;idim++)
      vdir[idim] *= (1.0/r);

    fwrite(&i,sizeof(type_int),1,pfrela);
    fwrite(&Seg[i].flag,sizeof(type_int),1,pfrela);    
    fwrite(&Gr[Seg[i].list[0]].Vnod[0],sizeof(type_real),3,pfrela);
    fwrite(&Gr[Seg[i].list[Seg[i].size-1]].Vnod[0],sizeof(type_real),3,pfrela);
    fwrite(&vdir[0],sizeof(type_real),3,pfrela);

  } //FINALIZA EL PARALELO

  fclose(pfrela);

  ///////////////////////////////////////////////////////

  BLUE("******************************\n");
}

#endif

extern void propiedades(const type_int NNN, const type_real *fof)
{
  char filename[200], name[200];
  type_int i, j, l, idim, Tid, *c;
  type_int totfil = 0;
  type_int binsper;
  type_real *rcil2;                 // En Kpc                    
  type_real rsep;                   // En Kpc
  type_real rlong;                  // En Kpc
  type_real rlong_2;                // En Kpc
  type_real RMAX, RMIN, LMAX;
  type_real r;
  FILE *pfdens;
  #ifdef CALCULA_MEDIA
  FILE *pfvel;
  #endif
  #ifdef SAVEPART
  FILE **pfpart;
  #endif
  #ifdef EXTEND
  FILE *pfextend;
  #endif
  struct dat_struct *dat;
  const type_real cut_len = 5000.;            // en Kpc

  dat = (struct dat_struct *) malloc(cp.nseg*sizeof(struct dat_struct));

  #ifdef EXTEND
    idim = nbins%2==0 ? 2*nbins : (type_int)(2.f*floor((float)nbins/2));
  #else
    idim = nbins;
  #endif

  l = 0;
  for(i=0;i<cp.nseg;i++)
  {
    if(Seg[i].flag != 2) continue;

    if((Seg[i].len < cut_len-1000.0) || (Seg[i].len > cut_len+1000.0)) continue;

    dat[l].nodo        = (type_real *)  malloc(idim*sizeof(type_real));
    dat[l].rbin_centre = (type_real *)  malloc(ncil *sizeof(type_real));

    dat[l].npart       = (type_int  **) malloc(ncil *sizeof(type_int  *));
    dat[l].rho_delta   = (type_real **) malloc(ncil *sizeof(type_real *));
    dat[l].mean        = (type_real **) malloc(ncil *sizeof(type_real *));
    dat[l].rms         = (type_real **) malloc(ncil *sizeof(type_real *));
    dat[l].mean_par    = (type_real **) malloc(ncil *sizeof(type_real *));
    dat[l].rms_par     = (type_real **) malloc(ncil *sizeof(type_real *));
    dat[l].mean_perp   = (type_real **) malloc(ncil *sizeof(type_real *));
    dat[l].rms_perp    = (type_real **) malloc(ncil *sizeof(type_real *));

    for(j=0;j<ncil;j++)
    {
      dat[l].npart[j]     = (type_int * ) malloc(idim*sizeof(type_int ));
      dat[l].rho_delta[j] = (type_real *) malloc(idim*sizeof(type_real));
      dat[l].mean[j]      = (type_real *) malloc(idim*sizeof(type_real));
      dat[l].rms[j]       = (type_real *) malloc(idim*sizeof(type_real));
      dat[l].mean_par[j]  = (type_real *) malloc(idim*sizeof(type_real));
      dat[l].rms_par[j]   = (type_real *) malloc(idim*sizeof(type_real));
      dat[l].mean_perp[j] = (type_real *) malloc(idim*sizeof(type_real));
      dat[l].rms_perp[j]  = (type_real *) malloc(idim*sizeof(type_real));
    }

    Seg[l] = Seg[i];
  
    l++;
  }

  dat = (struct dat_struct *) realloc(dat,l*sizeof(struct dat_struct));
  Seg = (struct segmentstd *) realloc(Seg,l*sizeof(struct segmentstd));
  cp.nseg = l;

  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(filename,"lbox %g Mpc %g Kpc\n",cp.lbox/1000.0f,cp.lbox);GREEN(filename);
  sprintf(filename,"cut LEN %f Mpc REALOCATEA %d fil\n",cut_len/1000.0f,cp.nseg);GREEN(filename);
  GREEN("**********************************\n");

  stadistic(cp.nseg,&RMAX,&RMIN,&LMAX); 

  #ifdef FIXED_SEPARATION
    ncil = (type_int)(RLEN/RSEP);
    RLEN *= 1000.f;
    RSEP *= 1000.f;
    assert(RLEN<RMAX);
    sprintf(message,"Num cylindres:          %d\n",ncil);BLUE(message);
  #endif

  fprintf(stdout,"RMAX %f Mpc\n",RMAX/1000.);
  fprintf(stdout,"RMIN %f Mpc\n",RMIN/1000.);
  fprintf(stdout,"LMAX %f Mpc\n",LMAX/1000.);
  binsper = ncil+1;

  //////////////////////////////////////////////////
  fprintf(stdout,"Build grid\n");

  grid.nobj = cp.npart;
 
  //if(LMAX>RMAX){
  //  grid.ngrid = (int)(cp.lbox/LMAX);
  //}else{ 
  //  grid.ngrid = (int)(cp.lbox/RMAX);
  //}

  grid.ngrid = (type_int)(cp.lbox/2000.);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }

  grid_init();
  grid_build();
  //////////////////////////////////////////////////

  c = (type_int *) malloc(NTHREADS*sizeof(type_int));
  #ifdef SAVEPART
    pfpart = (FILE **) malloc(NTHREADS*sizeof(FILE));
  #endif

  set_name(name,NNN,"densidad");
  sprintf(filename,"%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pfdens = fopen(filename,"w");

  #ifdef EXTEND
  set_name(name,NNN,"extend");
  sprintf(filename,"%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pfextend = fopen(filename,"w");
  #endif

  #ifdef CALCULA_MEDIA
  set_name(name,NNN,"vmedia");
  sprintf(filename,"../%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pfvel = fopen(filename,"w");
  #endif

  for(i=0;i<NTHREADS;i++)
  {
    c[i] = 0;
    #ifdef SAVEPART
      set_name(name,NNN,"particulas");
      sprintf(filename,"%.2d_%s_%.2f_%.2f.%.2d.bin",snap.num,name,fof[0],fof[1],i);
      pfpart[i] = fopen(filename,"w");
      fwrite(&c[i],sizeof(type_int),1,pfpart[i]);    
    #endif
  }

  #ifdef CALCULA_MEDIA
  
  BLUE("**** CALCULA VEL MEDIAS ******\n");
  fflush(stdout);

  #ifdef FIXED_SEPARATION

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(i, \
    rcil2,r,rsep,rlong,rlong_2) \
    shared(Gr,Seg,cp,ncil,binsper,c,RLEN,RMAX,RMIN,nbins)

  #else

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(i, \
    rcil2,r,rsep,rlong,rlong_2) \
    shared(Gr,Seg,cp,ncil,binsper,c,RMAX,RMIN,nbins)

  #endif
  for(i=0;i<cp.nseg;i++)
  {

    Seg[i].Vmedia = (type_real *) calloc(3*ncil,sizeof(type_real));
    rcil2  = (type_real *) malloc(binsper*sizeof(type_real));

    rsep    = Seg[i].len/(type_real)(nbins-1);
    rlong   = 2.0*rsep; 
    rlong_2 = rlong*0.5;

    #ifdef FIXED_SEPARATION
      r = RLEN; // RLEN<RMAX
    #else
      r = Seg[i].len;
    #endif

    if(r>RMAX)
    {
    #ifdef BIN_LOG
      logspace(rcil2,RMAX,RMIN,binsper);
    #else
      linspace(rcil2,RMAX,0.0,binsper);
    #endif
    }else{
    #ifdef BIN_LOG
      logspace(rcil2,r,r/50.,binsper);
    #else
      linspace(rcil2,r,0.0,binsper);
    #endif
    }

    cylinder_mean(i,rsep,rlong,rlong_2,rcil2);

    free(rcil2);

  } //FINALIZA EL PARALELO

  BLUE("******************************\n");

  set_name(name,NNN,"vmedia");
  sprintf(filename,"../%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pfvel = fopen(filename,"w");

  j = 0;
  fwrite(&j,sizeof(type_int),1,pfvel);        
  fwrite(&ncil,sizeof(type_int),1,pfvel);

  for(i=0;i<cp.nseg;i++)
  {

    if(Seg[i].flag != 2) continue;

    fwrite(&i,sizeof(int),1,pfvel);
    fwrite(&Seg[i].Vmedia[0],sizeof(type_real),3*ncil,pfvel);
    j++;

  }

  rewind(pfvel);
  fwrite(&j,sizeof(type_int),1,pfvel);
  fclose(pfvel);

  #endif

  BLUE("********** CALCULA ***********\n");
  fflush(stdout);

  #ifdef SAVEPART
  
    #ifdef FIXED_SEPARATION

      #pragma omp parallel for num_threads(NTHREADS) \
      schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
      r,rsep,rlong,rlong_2) shared(pfpart,Gr,P,Seg,cp, \
      RLEN,RMAX,RMIN,c,binsper,ncil,nbins,dat,stdout) reduction(+:totfil)

    #else

      #pragma omp parallel for num_threads(NTHREADS) \
      schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
      r,rsep,rlong,rlong_2) shared(pfpart,Gr,P,Seg,cp, \
      RMAX,RMIN,c,binsper,ncil,nbins,dat,stdout) reduction(+:totfil)

    #endif

  #else
  
   #ifdef FIXED_SEPARATION
   
     #pragma omp parallel for num_threads(NTHREADS) \
     schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
     r,rsep,rlong,rlong_2) shared(Gr,P,Seg,cp, \
     RLEN,RMAX,RMIN,c,binsper,ncil,nbins,dat,stdout) reduction(+:totfil)

   #else

     #pragma omp parallel for num_threads(NTHREADS) \
     schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
     r,rsep,rlong,rlong_2) shared(Gr,P,Seg,cp, \
     RMAX,RMIN,c,binsper,ncil,nbins,dat,stdout) reduction(+:totfil)
   
   #endif

  #endif
  for(i=0;i<cp.nseg;i++)
  {
    type_real  aux_mean, aux_rms;
    type_real *vmean_par, *vquad_par;
    type_real *vmean_perp, *vquad_perp;
    type_real **vmean, **vquad;
    struct node_sph **root;
    #ifdef SAVEPART
    type_int size_part = 0;
    type_int *vpart;
    vpart = (type_int *) calloc(cp.npart/32+1,sizeof(type_int));
    #endif

    rcil2  = (type_real *) malloc(binsper*sizeof(type_real));

    rsep    = Seg[i].len/(type_real)(nbins-1);
    rlong   = 2.0*rsep; 
    rlong_2 = rlong*0.5;

    #ifdef FIXED_SEPARATION
      r = RLEN; // RLEN<RMAX
    #else
      r = Seg[i].len;      
    #endif

    if(r>RMAX)
    {
      #ifdef BIN_LOG
        logspace(rcil2,RMAX,RMIN,binsper);
      #else
        linspace(rcil2,RMAX,0.0,binsper);
      #endif
    }else{
      #ifdef BIN_LOG
        logspace(rcil2,r,r/50.,binsper);
      #else
        linspace(rcil2,r,0.0,binsper);
      #endif
    }

    Tid = omp_get_thread_num();

    vmean_par = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_par = (type_real *) malloc(ncil*sizeof(type_real));
    vmean_perp = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_perp = (type_real *) malloc(ncil*sizeof(type_real));   

    root = (struct node_sph **) malloc(ncil*sizeof(struct node_sph *));
    vmean = (type_real **) malloc(ncil*sizeof(type_real *));
    vquad = (type_real **) malloc(ncil*sizeof(type_real *));
    for(j=0;j<ncil;j++)
    {
      root[j] = (struct node_sph *) malloc(nbins*sizeof(struct node_sph));
      vmean[j] = (type_real *) malloc(3*sizeof(type_real));
      vquad[j] = (type_real *) malloc(3*sizeof(type_real));
    }

    #ifdef SAVEPART
      calcular_mean(i,binsper,rsep,rlong,rlong_2,rcil2,\
      vmean,vquad,vmean_par,vquad_par,\
      vmean_perp,vquad_perp,dat[i].nodo,root,\
      &size_part,vpart);
    #else  
      calcular_mean(i,binsper,rsep,rlong,rlong_2,rcil2,\
      vmean,vquad,vmean_par,vquad_par,\
      vmean_perp,vquad_perp,dat[i].nodo,root);
    #endif  

    for(j=0;j<ncil;j++)
    {

      //r  = sqrt(rcil2[j]/rcil2[0]);
      //r += sqrt(rcil2[j+1]/rcil2[0]);
      //r *= 0.5;

      r  = sqrt(rcil2[j]);
      r += sqrt(rcil2[j+1]);
      r *= 0.5;

      //fwrite(&r,sizeof(type_real),1,pfdens[Tid]);       // rbin_centre
      dat[i].rbin_centre[j] = r;

      for(l=0;l<nbins;l++)
      {
        root[j][l].b -= 1.0;

        //fwrite(&root[j][l].a,sizeof(type_int),1,pfdens[Tid]);       // npart
        //fwrite(&root[j][l].b,sizeof(type_real),1,pfdens[Tid]); // pho - overdense
        dat[i].npart[j][l]     = root[j][l].a; // npart
        dat[i].rho_delta[j][l] = root[j][l].b; // rho - delta

        ////////////////////////////////////////////////////////////////////////////////////

        for(idim=0;idim<3;idim++)
        {
          root[j][l].c[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
          //root[j][l].d[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        }
        root[j][l].e *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        //root[j][l].f *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].g *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        //root[j][l].h *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;

        aux_mean = aux_rms = 0.0;

        for(idim=0;idim<3;idim++)
        {
          aux_mean  += root[j][l].c[idim];
          aux_rms   += root[j][l].d[idim] - (type_real)root[j][l].a*root[j][l].c[idim]*root[j][l].c[idim];
          //aux_rms   += root[j][l].d[idim] - root[j][l].c[idim]*root[j][l].c[idim];
        }
        
        aux_mean /= 3.0;        
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms)*(1.f/(type_real)(root[j][l].a-1)))/3.0 : 0.0;
        //aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms))/3.0 : 0.0;

        //fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean
        //fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms
        dat[i].mean[j][l] = aux_mean; // mean
        dat[i].rms[j][l]  = aux_rms;  // rms

        aux_mean = root[j][l].e;
        //aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].f - root[j][l].e*root[j][l].e)) : 0.0;
        aux_rms  = root[j][l].a>10 ? sqrt((fabs(root[j][l].f - \
        (type_real)root[j][l].a*root[j][l].e*root[j][l].e))*(1.f/(type_real)(root[j][l].a-1))) : 0.0;

        //fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean par
        //fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms  par
        dat[i].mean_par[j][l] = aux_mean; // mean par 
        dat[i].rms_par[j][l]  = aux_rms;  // rms  par

        aux_mean = root[j][l].g;
        //aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].h - root[j][l].g*root[j][l].g)) : 0.0;
        aux_rms  = root[j][l].a>10 ? sqrt((fabs(root[j][l].h - \
        (type_real)root[j][l].a*root[j][l].g*root[j][l].g))*(1.f/(type_real)(root[j][l].a-1))) : 0.0;

        //fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean perp        
        //fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms  perp
        dat[i].mean_perp[j][l] = aux_mean;// mean perp        
        dat[i].rms_perp[j][l]  = aux_rms; // rms  perp
      }

      free(vmean[j]);
      free(vquad[j]);
      free(root[j]);
    }

    free(root);
    free(vmean);
    free(vquad);
    free(vmean_par);
    free(vquad_par);
    free(vmean_perp);
    free(vquad_perp);

    #ifdef EXTEND

      #ifdef SAVEPART
        ext_mean(i,binsper,rsep,rlong,rlong_2,rcil2,&size_part,vpart,dat);
      #else
        ext_mean(i,binsper,rsep,rlong,rlong_2,rcil2,dat);
      #endif

    #endif

    #ifdef SAVEPART
    fwrite(&i,sizeof(type_int),1,pfpart[Tid]);
    fwrite(&size_part,sizeof(type_int),1,pfpart[Tid]);

    for(l=0;l<cp.npart;l++)
    {
      if(TestBit(vpart,l))    
        fwrite(&P[l].id,sizeof(type_int),1,pfpart[Tid]);
    }

    free(vpart);
    #endif

    c[Tid]++;
    totfil++;
    free(rcil2);

  } //FINALIZA EL PARALELO

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  

    #ifdef SAVEPART
      rewind(pfpart[i]);
      fwrite(&c[i],sizeof(type_int),1,pfpart[i]);
      fclose(pfpart[i]);
    #endif
  }

  fwrite(&cp.nseg,sizeof(type_int),1,pfdens);
  #ifdef EXTEND
    fwrite(&cp.nseg,sizeof(type_int),1,pfextend);
    idim = nbins%2==0 ? nbins/2 : (type_int)(floor((float)nbins/2));
  #endif

  for(i=0;i<cp.nseg;i++)
  {
    fwrite(&i,sizeof(int),1,pfdens);                   // Num Filamento
    fwrite(&Seg[i].flag,sizeof(int),1,pfdens);         // escribe la bandera
    fwrite(&Seg[i].len,sizeof(type_real),1,pfdens);    // escribe la longitud del segmento

    r = Gr[Seg[i].list[0]].NumPart*cp.Mpart;
    fwrite(&r,sizeof(type_real),1,pfdens);             // escribe la masa del primer nodo
    r = Gr[Seg[i].list[Seg[i].size-1]].NumPart*cp.Mpart;
    fwrite(&r,sizeof(type_real),1,pfdens);             // escribe la masa del ultimo nodo    

    fwrite(&nbins,sizeof(int),1,pfdens);               // la cantidad de cilindros que hizo

    for(j=0;j<nbins;j++)
     fwrite(&dat[i].nodo[j],sizeof(type_real),1,pfdens);

    for(j=0;j<ncil;j++)
    { 
      fwrite(&dat[i].rbin_centre[j],sizeof(type_real),1,pfdens);       // rbin_centre
      for(l=0;l<nbins;l++)
      {
          fwrite(&dat[i].npart[j][l],     sizeof(int),1,pfdens);       // npart
          fwrite(&dat[i].rho_delta[j][l], sizeof(type_real),1,pfdens); // rho - delta

          fwrite(&dat[i].mean[j][l], sizeof(type_real),1,pfdens); // mean
          fwrite(&dat[i].rms[j][l],  sizeof(type_real),1,pfdens); // rms

          fwrite(&dat[i].mean_par[j][l], sizeof(type_real),1,pfdens); // mean par
          fwrite(&dat[i].rms_par[j][l],  sizeof(type_real),1,pfdens); // rms  par

          fwrite(&dat[i].mean_perp[j][l], sizeof(type_real),1,pfdens); // mean par
          fwrite(&dat[i].rms_perp[j][l],  sizeof(type_real),1,pfdens); // rms  par
      }

      #ifndef EXTEND
        free(dat[i].npart[j]);
        free(dat[i].rho_delta[j]);
        free(dat[i].mean[j]);
        free(dat[i].rms[j]);
        free(dat[i].mean_par[j]);
        free(dat[i].rms_par[j]);
        free(dat[i].mean_perp[j]);
        free(dat[i].rms_perp[j]);
      #endif
    }

    #ifdef EXTEND
      j = 2*idim;
      fwrite(&i,sizeof(int),1,pfextend);                   // Num Filamento
      fwrite(&Seg[i].len,sizeof(type_real),1,pfextend);    // escribe la longitud del segmento
      fwrite(&j,sizeof(int),1,pfextend);                   // la cantidad de cilindros que voy a hacer

      /// INICIO ///
      for(j=nbins;j<nbins+idim;j++)
        fwrite(&dat[i].nodo[j],sizeof(type_real),1,pfextend);

      for(j=0;j<ncil;j++)
      { 
        fwrite(&dat[i].rbin_centre[j],sizeof(type_real),1,pfextend);       // rbin_centre
        for(l=nbins;l<nbins+idim;l++)
        {
          fwrite(&dat[i].npart[j][l],     sizeof(int),1,pfextend);       // npart
          fwrite(&dat[i].rho_delta[j][l], sizeof(type_real),1,pfextend); // rho - delta

          fwrite(&dat[i].mean[j][l], sizeof(type_real),1,pfextend); // mean
          fwrite(&dat[i].rms[j][l],  sizeof(type_real),1,pfextend); // rms

          fwrite(&dat[i].mean_par[j][l], sizeof(type_real),1,pfextend);  // mean par
          fwrite(&dat[i].rms_par[j][l],  sizeof(type_real),1,pfextend); // rms  par

          fwrite(&dat[i].mean_perp[j][l], sizeof(type_real),1,pfextend); // mean par
          fwrite(&dat[i].rms_perp[j][l],  sizeof(type_real),1,pfextend); // rms  par
        }
      }

      /// FINAL ///
      for(j=nbins+idim;j<nbins+2*idim;j++)
       fwrite(&dat[i].nodo[j],sizeof(type_real),1,pfextend);

      for(j=0;j<ncil;j++)
      { 
        fwrite(&dat[i].rbin_centre[j],sizeof(type_real),1,pfextend);     // rbin_centre

        for(l=nbins+idim;l<nbins+2*idim;l++)
        {
          fwrite(&dat[i].npart[j][l],     sizeof(int),1,pfextend);       // npart
          fwrite(&dat[i].rho_delta[j][l], sizeof(type_real),1,pfextend); // rho - delta

          fwrite(&dat[i].mean[j][l], sizeof(type_real),1,pfextend); // mean
          fwrite(&dat[i].rms[j][l],  sizeof(type_real),1,pfextend); // rms

          fwrite(&dat[i].mean_par[j][l], sizeof(type_real),1,pfextend);  // mean par
          fwrite(&dat[i].rms_par[j][l],  sizeof(type_real),1,pfextend); // rms  par

          fwrite(&dat[i].mean_perp[j][l], sizeof(type_real),1,pfextend); // mean par
          fwrite(&dat[i].rms_perp[j][l],  sizeof(type_real),1,pfextend); // rms  par
        }

        free(dat[i].npart[j]);
        free(dat[i].rho_delta[j]);
        free(dat[i].mean[j]);
        free(dat[i].rms[j]);
        free(dat[i].mean_par[j]);
        free(dat[i].rms_par[j]);
        free(dat[i].mean_perp[j]);
        free(dat[i].rms_perp[j]);
      }
    #endif

    free(dat[i].nodo);
    free(dat[i].rbin_centre);
    free(dat[i].npart);
    free(dat[i].rho_delta);
    free(dat[i].mean);
    free(dat[i].rms);
    free(dat[i].mean_par);
    free(dat[i].rms_par);
    free(dat[i].mean_perp);
    free(dat[i].rms_perp);

  } 

  fclose(pfdens);
  #ifdef EXTEND
  fclose(pfextend);
  #endif

  free(dat);
  free(c);
  grid_free();

  BLUE("******************************\n");
  fprintf(stdout,"Calcula %d Filamentos\n",totfil);  
  BLUE("******************************\n");

  return;
}
