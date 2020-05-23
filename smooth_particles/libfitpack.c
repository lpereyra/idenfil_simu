#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "colores.h"
#include "libfitpack.h"

#define IX(i,j) 3*i+j

static void smooth(int n, int idim, float *Pos, double s, int npoint)
{
  int k = 3, ipar = 0, iopt = 0, m = 0;  
  int i,j,ier, nest, nx, nc, lwrk;
  int *iwrk;
  double ub, ue, fp;
  double *c, *t, *wrk, *wa;
  double *x, *w, *u;

  assert(n > k);
  assert(k>=1 && k<=5);

  nest = n+2*k>2*k+3 ? n+2*k : 2*k+3;

  nx = idim*n;
  nc = idim*nest;
  lwrk = n*(k + 1) + nest*(6 + idim + 3*k);
  j = nc + 2*nest + lwrk;

  wa = (double *) malloc(j*sizeof(double));
  t = wa;
  c = t + nest;
  wrk = c + nc;
  iwrk = (int *)(wrk + lwrk);

  x = (double *) malloc(nx*sizeof(double));
  w = (double *) malloc(n*sizeof(double));
  u = (double *) malloc(n*sizeof(double));
  
  for(j=0;j<n;j++)
  {
    for(i=0;i<idim;i++)
    {
      x[idim*j+i] = (double)Pos[idim*j+i];
    }
    w[j] = 1.0f;
    u[j] = 0.0f;
  }

  parcur_(&iopt, &ipar, &idim, &n, u, &nx, x, w, &ub, &ue, &k, \
          &s, &nest, &m, t, &nc, c, &fp, wrk, &lwrk, iwrk, &ier);

  if(ier == 10) 
  {
    assert(0);
  }

  if(ier > 0 && m == 0) {
     n = 1;
  }

  free(u);
  free(x);
  free(w);
  
  nx = idim*npoint;
  x   = (double *) malloc(npoint*sizeof(double));
  w   = (double *) malloc(npoint*sizeof(double));

  fp = 1./(npoint-1);
  for(j=0;j<npoint;j++)
    x[j] = j*fp;

  for(i=0;i<idim;i++)
  {
    splev_(t,&m,c+i*m,&k,x,w,&npoint,&iopt,&ier);

    for(j=0;j<npoint;j++)
      Pos[j*idim+i] = (float)w[j];

    if(ier == 10) 
      assert(0);
  }

  free(x);
  free(w);
  free(wa);

  return;
}

extern void suaviza()
{
  type_int i,k,idim;
  type_real diff, *Pos_aux;
  #ifdef ITERA
  type_int j;
  #endif
  type_real fix_array[6];

  cp.lbox /= POSFACTOR;

  for(i=0;i<cp.nseg;i++)
    for(k=0;k<Seg[i].size;k++)
      for(idim=0;idim<3;idim++)
        Seg[i].Pos_list[IX(k,idim)] /= POSFACTOR;

  #ifdef ITERA

  GREEN("*********** ITERA ************\n");
  GREEN("******************************\n");
  fflush(stdout); 

  for(i=0;i<cp.nseg;i++)
  {
    if(Seg[i].size==2) continue;
      
    Pos_aux = (type_real *) malloc(3*Seg[i].size*sizeof(type_real));

    for(j=0;j<NITERA;j++)
    {
      for(k=0;k<Seg[i].size;k++)
      {      
        for(idim=0;idim<3;idim++)
        {
          Pos_aux[IX(k,idim)] = Seg[i].Pos_list[IX(k,idim)];

          diff = Seg[i].Pos_list[IX(k,idim)] - Seg[i].Pos_list[IX(0,idim)];
          #ifdef PERIODIC
          if(diff> 0.5*cp.lbox) Pos_aux[IX(k,idim)] -= cp.lbox;
          if(diff<-0.5*cp.lbox) Pos_aux[IX(k,idim)] += cp.lbox;
          #endif
        }
      }

      for(k=1;k<Seg[i].size-1;k++)
      {      
        for(idim=0;idim<3;idim++)
        {
          Seg[i].Pos_list[IX(k,idim)] = \
          Pos_aux[IX(k,idim)];

          //Seg[i].Pos_list[IX(k,idim)] = \
          //0.25f*Pos_aux[IX((k-1),idim)]  + \
          //0.50f*Pos_aux[IX(k,idim)]      + \
          //0.25f*Pos_aux[IX((k+1),idim)];
  
          #ifdef PERIODIC
          if(Seg[i].Pos_list[IX(k,idim)] >= cp.lbox) 
            Seg[i].Pos_list[IX(k,idim)] -= cp.lbox;

          if(Seg[i].Pos_list[IX(k,idim)] < 0.0f) 
            Seg[i].Pos_list[IX(k,idim)] += cp.lbox;
          #endif
        }       
      }

    }

    free(Pos_aux);
  }

  GREEN("********* END ITERA **********\n");
  GREEN("******************************\n");
  fflush(stdout); 
  
  #endif

  GREEN("********** SUAVIZA ***********\n");
  GREEN("******************************\n");
  fflush(stdout); 

  for(i=0;i<cp.nseg;i++)
  {
    #ifdef FIX_NSMOOTH
      type_int Nsm = NSMOOTH;
    #else
      type_int Nsm = ceil(Seg[i].len/(type_real)RSPACE);
    #endif

    if(Seg[i].size>=4)
    {

      if(Nsm<=Seg[i].size) Nsm = Seg[i].size;

      Seg[i].Pos_list = (type_real *) realloc(Seg[i].Pos_list,3*Nsm*sizeof(type_real));

      for(k=0;k<Seg[i].size;k++)
      {      
        for(idim=0;idim<3;idim++)
        {
          diff = Seg[i].Pos_list[IX(k,idim)] - Seg[i].Pos_list[IX(0,idim)];
          #ifdef PERIODIC
          if(diff> 0.5*cp.lbox) Seg[i].Pos_list[IX(k,idim)] -= cp.lbox;
          if(diff<-0.5*cp.lbox) Seg[i].Pos_list[IX(k,idim)] += cp.lbox;
          #endif
        }
      }

      for(idim=0;idim<3;idim++)
      {
        fix_array[idim]   = Seg[i].Pos_list[IX(0,idim)];
        fix_array[idim+3] = Seg[i].Pos_list[IX((Seg[i].size-1),idim)];
      }

      smooth((int)Seg[i].size,3,Seg[i].Pos_list,2.0,(int)Nsm);

      Seg[i].size = Nsm;

      for(idim=0;idim<3;idim++)
      {
        Seg[i].Pos_list[IX(0,idim)]               = fix_array[idim];
        Seg[i].Pos_list[IX((Seg[i].size-1),idim)] = fix_array[idim+3];
      }
     
      for(k=0;k<Seg[i].size;k++)
      {      
        for(idim=0;idim<3;idim++)
        { 
          #ifdef PERIODIC
          Seg[i].Pos_list[IX(k,idim)] = Seg[i].Pos_list[IX(k,idim)] >= cp.lbox ? \
            Seg[i].Pos_list[IX(k,idim)]-cp.lbox : Seg[i].Pos_list[IX(k,idim)];

          Seg[i].Pos_list[IX(k,idim)] = Seg[i].Pos_list[IX(k,idim)] < 0.0f ? \
            Seg[i].Pos_list[IX(k,idim)]+cp.lbox : Seg[i].Pos_list[IX(k,idim)];
          #endif
        }
      }
    }

    for(k=0;k<Seg[i].size;k++)
      for(idim=0;idim<3;idim++)
        Seg[i].Pos_list[IX(k,idim)] *= POSFACTOR;
  }

  cp.lbox *= POSFACTOR;

  GREEN("******** END SUAVIZADO ********\n");
  GREEN("*******************************\n");
  fflush(stdout); 

  GREEN("***** CALCULA PROPIEDADES *****\n");
  GREEN("*******************************\n");
  fflush(stdout); 

  #pragma omp parallel for num_threads(NTHREADS) default(none) \
  private(i,k,idim,diff,fix_array) shared(cp,Seg,stdout)
  for(i=0;i<cp.nseg;i++)
  {
    /////////////////////////////////////////////////////////////////////  

    Seg[i].len = 0.0f;
    for(k=1;k<Seg[i].size;k++)
    {
      diff = 0.0f;
      for(idim=0;idim<3;idim++)
      { 
        fix_array[idim] = Seg[i].Pos_list[IX(k,idim)]-Seg[i].Pos_list[IX((k-1),idim)];
 
        #ifdef PERIODIC
        if(fix_array[idim] >  0.5*cp.lbox) fix_array[idim] -= cp.lbox;
        if(fix_array[idim] < -0.5*cp.lbox) fix_array[idim] += cp.lbox;
        #endif

        diff += fix_array[idim]*fix_array[idim];
      }
      Seg[i].len += sqrt(diff);
    }

    /////////////////////////////////////////////////////////////////////  
    
    diff = 0.0f;
    for(idim=0;idim<3;idim++)
    { 
      fix_array[idim+3] = Seg[i].Pos_list[IX((Seg[i].size-1),idim)]-Seg[i].Pos_list[IX(0,idim)];
   
      #ifdef PERIODIC
      if(fix_array[idim+3] >  0.5*cp.lbox) fix_array[idim+3] -= cp.lbox;
      if(fix_array[idim+3] < -0.5*cp.lbox) fix_array[idim+3] += cp.lbox;
      #endif
  
      diff += fix_array[idim+3]*fix_array[idim+3];
    }

    diff = sqrt(diff);

    for(idim=0;idim<3;idim++)
      fix_array[idim+3] /= diff;
  
    Seg[i].elong = diff/Seg[i].len;

    ///////////////////////////////////////////////////////////////////// 
  
    /////////////////////////////////////////////////////////////////////  
    
    Seg[i].rms = 0.0f;
    for(k=1;k<Seg[i].size-1;k++)
    {
      diff = 0.0f;
      for(idim=0;idim<3;idim++)
      { 
        fix_array[idim] = Seg[i].Pos_list[IX(k,idim)]-Seg[i].Pos_list[IX(0,idim)];
 
        #ifdef PERIODIC
        if(fix_array[idim] >  0.5*cp.lbox) fix_array[idim] -= cp.lbox;
        if(fix_array[idim] < -0.5*cp.lbox) fix_array[idim] += cp.lbox;
        #endif

      }

      diff = fix_array[0]*fix_array[3]+fix_array[1]*fix_array[4]+fix_array[2]*fix_array[5];              // dot
      diff = fix_array[0]*fix_array[0]+fix_array[1]*fix_array[1]+fix_array[2]*fix_array[2] - diff*diff;  // dist - dot
      diff = fabs(diff);  // por errores de redondeo
      Seg[i].rms += diff;

    }

    Seg[i].rms /= (type_int)Seg[i].size;
    Seg[i].rms = sqrt(Seg[i].rms);

  }

  return;
}

