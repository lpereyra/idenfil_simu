#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "calcula.h"
#include "colores.h"
#include "leesnap.h"
#include "grid.h"
#include "list.h"

#define IX(i,j,n) n*j+i

static const type_real r_mean = 1500.; // En Kpc

static type_real **matrix_mean_per;
static type_real **matrix_quad_per;
static type_real **matrix_mean_par;
static type_real **matrix_quad_par;
static type_int  **matrix_npart;
static type_int  *matrix_middle;
#ifdef LOCK
  static omp_lock_t *lock; 
#endif
#ifdef CALCULA_MEDIA
  static type_int **matrix_sum_npart;
#endif

static void set_name(char * name, const type_int NNN, const char * prefix)
{
  #ifdef CUT_ELONGACION
    sprintf(name,"../test_elongacion/");
  #else
    sprintf(name,"./");
  #endif

  #ifdef ORIGINAL
    sprintf(name,"%s%.2d_%.4d",name,snap.num,NNN);
  #else
    sprintf(name,"%s%.2d_%.4d_smooth",name,snap.num,NNN);
  #endif

  #ifdef CUT_ELONGACION
    sprintf(name,"%s_%.2f",name,CUT_ELONG);
  #endif

  #ifdef NEW
    sprintf(name,"%s_new",name);
  #endif

  #ifdef MCRITIC
    sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%s",name,prefix);

  #ifdef CALCULA_VCM
    sprintf(name,"%s_vcm",name);
  #endif

  #ifdef CALCULA_MEDIA
    sprintf(name,"%s_media",name);
  #endif

  #ifdef EXTEND
    sprintf(name,"%s_extend",name);
  #endif

  #ifdef FIXED_SEPARATION
    sprintf(name,"%s_%.2f_%.2f_Mpc",name,RLEN,RSEP);
  #else
    sprintf(name,"%s_%.2d_%.2d_%.2f_Mpc",name,nbins,ncil,RAUX/1000.);
  #endif

  #ifdef BIN_LOG
    sprintf(name,"%s_LOG",name);
  #else
    sprintf(name,"%s_LIN",name);
  #endif

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);

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

static void calc_fil_dis(const type_int idpar, const type_int idfil, const type_real * restrict rcil2)
{
  int j, k, idim, ix, iy, iseg;
  type_real racum_min,rx_min,ry_min;
  type_real r,rsep,rlen,racum,dot,dis;
  type_real Pos_cent[3], vdir[3], vprima[3], delta[3];

  iseg = -1;
  dis = racum = 0.0f;
  ry_min = cp.lbox*cp.lbox;
  racum_min = rx_min = -cp.lbox;
  #ifdef EXTEND
  rlen = Seg[idfil].len_extend;
  #else
  rlen = Seg[idfil].len;
  #endif

  ////////////////////////////////////////////////////////////////////////////////////////////

  //
  //REGION EXTEND
  //
  
  rsep = RAUX;

  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    vdir[idim] = Gr[Seg[idfil].list[1]].Pos[idim]-Gr[Seg[idfil].list[0]].Pos[idim];    
    #ifdef PERIODIC
    if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
    if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
    #endif
    r += vdir[idim]*vdir[idim];
  }
  r = sqrt(r);

  for(idim=0;idim<3;idim++)
  {
    vdir[idim] *= (1/r); // directo
    Pos_cent[idim] = Gr[Seg[idfil].list[0]].Pos[idim] - rsep*vdir[idim]; // cambio el origen
    #ifdef PERIODIC
    Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
    Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
    #endif
  }

  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
    #ifdef PERIODIC
    if(delta[idim]> 0.5*cp.lbox) delta[idim] -= cp.lbox;
    if(delta[idim]<-0.5*cp.lbox) delta[idim] += cp.lbox;
    #endif
    r += delta[idim]*delta[idim];
  }

  dot = delta[0]*vdir[0] + delta[1]*vdir[1] + delta[2] * vdir[2];

  // SI ENTRA EN MI CILINDRO
  if(dot>=0.0f && dot<=rsep) 
  {
    dis = fabs(r - dot*dot);
    if(dis<ry_min && dis<rcil2[0])
    {
      ry_min = dis;
      rx_min = dot; 
      racum_min = racum; 
      iseg = 0;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  //
  // REGION INTERNA
  //

  for(k=1;k<Seg[idfil].size;k++)
  {

    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[idfil].list[k]].Pos[idim]-Gr[Seg[idfil].list[k-1]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);
    rsep = r; // para mirar si caigo adentro del cilindro

    for(idim=0;idim<3;idim++)
    {
      vdir[idim] *= (1/r); // directo
      Pos_cent[idim] = Gr[Seg[idfil].list[k-1]].Pos[idim];
      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
      #ifdef PERIODIC
      if(delta[idim]> 0.5*cp.lbox) delta[idim] -= cp.lbox;
      if(delta[idim]<-0.5*cp.lbox) delta[idim] += cp.lbox;
      #endif
      r += delta[idim]*delta[idim];
    }

    dot = delta[0]*vdir[0] + delta[1]*vdir[1] + delta[2] * vdir[2];

    // SI ENTRA EN MI CILINDRO
    if(dot>=0.0f && dot<=rsep) 
    {
      dis = fabs(r - dot*dot);
      if(dis<ry_min && dis<rcil2[0])
      {
        ry_min = dis;
        rx_min = dot; 
        racum_min = racum; 
        iseg = k;
      }
    }

    racum += rsep;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  
  //
  //REGION EXTEND
  //

  rsep   = RAUX;

  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    vdir[idim] = Gr[Seg[idfil].list[Seg[idfil].size-1]].Pos[idim]-Gr[Seg[idfil].list[Seg[idfil].size-2]].Pos[idim];    
    #ifdef PERIODIC
    if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
    if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
    #endif
    r += vdir[idim]*vdir[idim];
  }
  r = sqrt(r);

  for(idim=0;idim<3;idim++)
  {
    vdir[idim] *= (1/r); // directo
    Pos_cent[idim] = Gr[Seg[idfil].list[Seg[idfil].size-1]].Pos[idim];
    #ifdef PERIODIC
    Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
    Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
    #endif
  }

  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
    #ifdef PERIODIC
    if(delta[idim]> 0.5*cp.lbox) delta[idim] -= cp.lbox;
    if(delta[idim]<-0.5*cp.lbox) delta[idim] += cp.lbox;
    #endif
    r += delta[idim]*delta[idim];
  }

  dot = delta[0]*vdir[0] + delta[1]*vdir[1] + delta[2]*vdir[2];

  // SI ENTRA EN MI CILINDRO
  if(dot>=0.0f && dot<=rsep) 
  {
    dis = fabs(r - dot*dot);
    if(dis<ry_min && dis<rcil2[0])
    {
      ry_min = dis;
      rx_min = dot;
      racum_min = racum; 
      iseg = Seg[idfil].size;
    }
  }

  if(iseg<0)
    return;

  j = nbins; 

  #ifdef EXTEND
    j *= 2;
  #endif

  if(iseg==0)
  {
    return;

  }else if(iseg==Seg[idfil].size){

    return;

  }else{

    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[idfil].list[iseg]].Pos[idim]-Gr[Seg[idfil].list[iseg-1]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);

    for(idim=0;idim<3;idim++)
    {
      vdir[idim] *= (1/r); // directo
      Pos_cent[idim] = Gr[Seg[idfil].list[iseg-1]].Pos[idim];

      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

  }

  r = rlen/(type_real)j;
  rsep  = racum_min + rx_min;
  ix = (int)(rsep/r); // bin round
  if(ix>=j)
    return;

  iy = 0;
  for(k=1;k<ncil;k++)
    iy = ry_min<rcil2[k] ? k : iy;

  k = j;
  j = IX(ix,iy,j);
  assert(j<ncil*k);

  for(idim=0;idim<3;idim++)
  {
   vprima[idim] = P[idpar].Vel[idim]; // asigna

   #ifdef CALCULA_MEDIA
     vprima[idim] -= Seg[idfil].Vmedia[3*iy+idim]; // asigna
   #endif

   #ifdef CALCULA_VCM
     vprima[idim] -= Seg[idfil].Vmedia[idim]; // asigna
   #endif

   delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
   #ifdef PERIODIC
   if(delta[idim] >  0.5*cp.lbox) delta[idim] -= cp.lbox;
   if(delta[idim] < -0.5*cp.lbox) delta[idim] += cp.lbox;
   #endif
  }

  r = dot = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    dot += (vprima[idim]*vdir[idim]);
  
    delta[idim] -= rx_min*vdir[idim];
  
    #ifdef PERIODIC
    if(delta[idim]>  0.5f*cp.lbox) delta[idim] -= cp.lbox;
    if(delta[idim]< -0.5f*cp.lbox) delta[idim] += cp.lbox;
    #endif
  
    r += (delta[idim]*delta[idim]);
  }
  
  r = sqrt(r);

  racum = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    delta[idim] *= (1.0f/r);
    racum += (vprima[idim]*delta[idim]);
  }

  #ifdef EXTEND
    r = rsep-RAUX;
    r *= (1.0/Seg[idfil].len_extend);
  #else
    r = rsep;
    r *= (1.0/Seg[idfil].len);
  #endif


  if(r>0.25 && r<0.75)
  {
    if(r_mean<sqrt(ry_min))
      iseg = -1;
    else
      iseg =  1;
  }else{
    iseg = -1;
  }

  #ifdef LOCK
  omp_set_lock(&lock[idfil]);
  #endif
  matrix_mean_par[idfil][j] += dot;
  matrix_quad_par[idfil][j] += (dot*dot);
  matrix_mean_per[idfil][j] += racum;
  matrix_quad_per[idfil][j] += (racum*racum);  
  matrix_npart[idfil][j]++;
  if(iseg>0) 
    matrix_middle[idfil]++;
  #ifdef LOCK
  omp_unset_lock(&lock[idfil]);
  #endif

  return;
}

static void calc_search_fil(const type_real * restrict Pos_cent,  const type_real r_2, struct list ** lista)
{	
	type_int  i;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  type_real fac,dis;
  type_real Posprima[3];

  fac   = (type_real)grid.ngrid/cp.lbox;
  ixc  = (long)(Pos_cent[0]*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (long)(Pos_cent[1]*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (long)(Pos_cent[2]*fac);
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
          Posprima[0] = Gr[i].Pos[0] - Pos_cent[0];
          Posprima[1] = Gr[i].Pos[1] - Pos_cent[1];
          Posprima[2] = Gr[i].Pos[2] - Pos_cent[2];

          #ifdef PERIODIC
          if(Posprima[0] >  0.5f*cp.lbox) Posprima[0] = Posprima[0] - cp.lbox;
          if(Posprima[1] >  0.5f*cp.lbox) Posprima[1] = Posprima[1] - cp.lbox;
          if(Posprima[2] >  0.5f*cp.lbox) Posprima[2] = Posprima[2] - cp.lbox;
          if(Posprima[0] < -0.5f*cp.lbox) Posprima[0] = Posprima[0] + cp.lbox;
          if(Posprima[1] < -0.5f*cp.lbox) Posprima[1] = Posprima[1] + cp.lbox;
          if(Posprima[2] < -0.5f*cp.lbox) Posprima[2] = Posprima[2] + cp.lbox;
          #endif

          dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2];

          if(dis<r_2)
          {
            struct data_list dat;
            dat.id = Gr[i].save;
            dat.r  = dis;  // el cuadrado
            push_list(lista, dat);
          }// cierra el if

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

#ifdef CALCULA_MEDIA

  static void calc_fil_mean(const type_int idpar, const type_int idfil, const type_real * restrict rcil2)
  {
    int k, idim, iy, iseg;
    type_real r,rsep,racum,dot,dis,ry_min;
    type_real Pos_cent[3], vdir[3], delta[3];
    
    iseg = -1;
    dis = racum = 0.0f;
    ry_min = cp.lbox*cp.lbox;
    
    ////////////////////////////////////////////////////////////////////////////////////////////
  
    //
    // REGION INTERNA
    //
  
    for(k=1;k<Seg[idfil].size;k++)
    {
  
      r = 0.0f;
      for(idim=0;idim<3;idim++)
      {
        vdir[idim] = Gr[Seg[idfil].list[k]].Pos[idim]-Gr[Seg[idfil].list[k-1]].Pos[idim];    
        #ifdef PERIODIC
        if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
        if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
        #endif
        r += vdir[idim]*vdir[idim];
      }
      r = sqrt(r);
      rsep = r; // para mirar si caigo adentro del cilindro
  
      //#pragma unroll
      for(idim=0;idim<3;idim++)
      {
        vdir[idim] *= (1/r); // directo
        Pos_cent[idim] = Gr[Seg[idfil].list[k-1]].Pos[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }
  
      r = 0.0f;
      for(idim=0;idim<3;idim++)
      {
        delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
        #ifdef PERIODIC
        if(delta[idim]> 0.5*cp.lbox) delta[idim] -= cp.lbox;
        if(delta[idim]<-0.5*cp.lbox) delta[idim] += cp.lbox;
        #endif
        r += delta[idim]*delta[idim];
      }
  
      dot = delta[0]*vdir[0] + delta[1]*vdir[1] + delta[2] * vdir[2];
  
      // SI ENTRA EN MI CILINDRO
      if(dot>=0.0f && dot<=rsep) 
      {
        dis = fabs(r - dot*dot);
        if(dis<ry_min && dis<rcil2[0])
        {
          ry_min = dis;
          iseg = k;
        }
      }
  
      racum += rsep;
    }
  
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    if(iseg<0)
      return;
  
    iy = 0;
    for(k=1;k<ncil;k++)
      iy = ry_min<rcil2[k] ? k : iy;
  
    assert(iy<ncil);
    
    #ifdef LOCK
    omp_set_lock(&lock[idfil]);
    #endif
    Seg[idfil].Vmedia[3*iy + 0] += P[idpar].Vel[0];
    Seg[idfil].Vmedia[3*iy + 1] += P[idpar].Vel[1];
    Seg[idfil].Vmedia[3*iy + 2] += P[idpar].Vel[2];
    matrix_sum_npart[idfil][iy]++;
    #ifdef LOCK
    omp_unset_lock(&lock[idfil]);
    #endif
    
    return;
  }
  
  static void calcular_media(const type_int i, const type_real rbus)
  {
    struct list *root;
    struct list *aux_lista_fil = NULL;
    type_real *rcil2;
  
    calc_search_fil(P[i].Pos,rbus*rbus,&aux_lista_fil);             
    remove_duplicates(aux_lista_fil);
  
    rcil2  = (type_real *) malloc((ncil+1)*sizeof(type_real));
  
    #ifdef BIN_LOG
      logspace(rcil2,RAUX,R_INIT,ncil+1);
    #else
      linspace(rcil2,RAUX,R_INIT,ncil+1);
    #endif

    while(aux_lista_fil != NULL)
    {
      root = aux_lista_fil;
      calc_fil_mean(i,root->data.id,rcil2);
      aux_lista_fil = aux_lista_fil->next;
      free(root);
    }
  
    free(rcil2);
  
    return;
  }


#endif

static void calcular_pert(const type_int i, const type_real rbus)
{
  struct list *root;
  struct list *aux_lista_fil = NULL;
  type_real *rcil2;

  calc_search_fil(P[i].Pos,rbus*rbus,&aux_lista_fil);             

  remove_duplicates(aux_lista_fil);

  rcil2  = (type_real *) malloc((ncil+1)*sizeof(type_real));

  #ifdef BIN_LOG
    logspace(rcil2,RAUX,R_INIT,ncil+1);
  #else
    linspace(rcil2,RAUX,R_INIT,ncil+1);
  #endif

  while(aux_lista_fil != NULL)
  {
    root = aux_lista_fil;
    calc_fil_dis(i,root->data.id,rcil2);
    aux_lista_fil = aux_lista_fil->next;
    free(root);
  }

  free(rcil2);

  return;
}

extern void propiedades(type_int NNN, type_real *fof)
{
  type_int i,j,k;
  const type_real r_aux = 2.0f*RAUX; // RADIO DE BUSQUEDA EN Kpc
  type_real *rcil2, *volinv;
  type_real aux,r,rlong;
  char filename[200];
  FILE *pf;

  #ifdef EXTEND

    for(i=0;i<cp.nseg;i++)
    {
      Seg[i].len_extend = 0.0f;
      for(k=1;k<Seg[i].size;k++)
      {
        r = 0.0f;
        for(j=0;j<3;j++)
        {
          aux = Gr[Seg[i].list[k]].Pos[j]-Gr[Seg[i].list[k-1]].Pos[j];    
          #ifdef PERIODIC
          if(aux> 0.5*cp.lbox) aux -= cp.lbox;
          if(aux<-0.5*cp.lbox) aux += cp.lbox;
          #endif
          r += aux*aux;
        }
        Seg[i].len_extend += sqrt(r);
      }
    }

  #endif

  grid.nobj = cp.ngrup;
  grid.ngrid = (long)(cp.lbox/r_aux);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }else{
    fprintf(stdout,"Using NGRID = %lu\n",grid.ngrid);
  }

  fflush(stdout);

  grid_init();
  grid_build();

  j = nbins;
  #ifdef EXTEND
  j += nbins;
  #endif

  matrix_mean_per = (type_real **) malloc(cp.nseg*sizeof(type_real *));
  matrix_quad_per = (type_real **) malloc(cp.nseg*sizeof(type_real *));
  matrix_mean_par = (type_real **) malloc(cp.nseg*sizeof(type_real *));
  matrix_quad_par = (type_real **) malloc(cp.nseg*sizeof(type_real *));
  matrix_npart    = (type_int  **) malloc(cp.nseg*sizeof(type_int  *));
  matrix_middle   = (type_int  *)  malloc(cp.nseg*sizeof(type_int));
  #ifdef LOCK
    lock = (omp_lock_t *) malloc(cp.nseg*sizeof(omp_lock_t));
  #endif
  #ifdef CALCULA_MEDIA
    matrix_sum_npart = (type_int **) malloc(cp.nseg*sizeof(type_int *));
  #endif

  for(i=0;i<cp.nseg;i++)  
  {
    matrix_mean_per[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_quad_per[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_mean_par[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_quad_par[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_npart[i]    = (type_int  *) calloc(ncil*j,sizeof(type_int ));
    #ifdef LOCK
      omp_init_lock(&(lock[i]));
    #endif
    #ifdef CALCULA_MEDIA
      matrix_sum_npart[i] = (type_int  *) calloc(ncil  ,sizeof(type_int ));
      Seg[i].Vmedia       = (type_real *) calloc(3*ncil,sizeof(type_real));
    #endif  
    matrix_middle[i] = 0;
  }

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  BLUE("******************************\n");
  GREEN("******************************\n");
  fflush(stdout);

  #ifdef CALCULA_VCM

    GREEN("CALCULA VCM\n");

    for(i=0;i<cp.nseg;i++)  
    {
      for(k=0;k<3;k++)
      {
        Seg[i].Vmedia[k] = (Seg[i].Mass[0]*Seg[i].Vnodos[k]+    \
                            Seg[i].Mass[1]*Seg[i].Vnodos[k+3])/ \
                           (Seg[i].Mass[0]+Seg[i].Mass[1]);
      }
    }

    GREEN("Termina el calculo\n");

  #endif

  #ifdef CALCULA_MEDIA

    GREEN("CALCULA MEDIA\n");

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) \
    private(i) shared(cp,stdout)
    for(i=0;i<cp.npart;i++)
    {
      if(i%1000000==0) fprintf(stdout,"%u %.4f\n",i,(float)i/(float)cp.npart);
      calcular_media(i,r_aux);
    }

    for(i=0;i<cp.nseg;i++)  
    {
      for(k=0;k<ncil;k++)
      {
        Seg[i].Vmedia[3*k + 0] *= (matrix_sum_npart[i][k]>0 ? 1.0f/(type_real)matrix_sum_npart[i][k] : 0.0);
        Seg[i].Vmedia[3*k + 1] *= (matrix_sum_npart[i][k]>0 ? 1.0f/(type_real)matrix_sum_npart[i][k] : 0.0);
        Seg[i].Vmedia[3*k + 2] *= (matrix_sum_npart[i][k]>0 ? 1.0f/(type_real)matrix_sum_npart[i][k] : 0.0);
      }

      free(matrix_sum_npart[i]);
    }

    free(matrix_sum_npart);

    GREEN("Termina el calculo\n");

  #endif

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) \
  private(i) shared(cp,stdout)
  for(i=0;i<cp.npart;i++)
  {
    if(i%1000000==0) fprintf(stdout,"%u %.4f\n",i,(float)i/(float)cp.npart);
    calcular_pert(i,r_aux);
  }

  GREEN("******************************\n");
  BLUE("******************************\n");
  fflush(stdout);

  grid_free();

  rcil2  = (type_real *) malloc((ncil+1)*sizeof(type_real));
  volinv = (type_real *) malloc(ncil*sizeof(type_real));

  set_name(filename,NNN,"matrix");
  pf = fopen(filename,"w");

  fwrite(&cp.nseg,sizeof(type_int),1,pf);        
  fwrite(&j,sizeof(type_int),1,pf);        
  fwrite(&ncil,sizeof(type_int),1,pf);        

  for(i=0;i<cp.nseg;i++) 
  {
    fwrite(&i,sizeof(type_int),1,pf);                   // Num Filamento
    fwrite(&Seg[i].flag,sizeof(type_int),1,pf);         // escribe la bandera
    fwrite(&Seg[i].len,sizeof(type_real),1,pf);    // escribe la longitud del segmento

    fwrite(&Seg[i].Mass[0],sizeof(type_real),1,pf);             // escribe la masa del primer nodo
    fwrite(&Seg[i].Mass[1],sizeof(type_real),1,pf);             // escribe la masa del ultimo nodo   

    aux  = (type_real)matrix_middle[i]*(cp.lbox*cp.lbox*cp.lbox);
    aux /= (type_real)cp.npart;
    aux /= (Seg[i].len*0.5*4.0*M_PI*r_mean*r_mean);
    aux -= 1.0f;
    fwrite(&aux,sizeof(type_real),1,pf);             // escribe la masa dentro de un cilindro de r_mean

    #ifdef BIN_LOG
      logspace(rcil2,RAUX,R_INIT,ncil+1);
    #else
      linspace(rcil2,RAUX,R_INIT,ncil+1);
    #endif

    #ifdef EXTEND
    rlong = Seg[i].len_extend;
    #else
    rlong = Seg[i].len;
    #endif

    rlong /= (type_real)j;

    for(k=0;k<ncil;k++)
    {
      volinv[k] = M_PI*rlong*rcil2[k];
    
      #ifdef HUECOS
        r = M_PI*rlong*rcil2[k+1];
        if(j != ncil-1) // para hacer el ultimo cilindro entero
          volinv[k] -= r;
      #endif  

      volinv[k] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[k])/(type_real)cp.npart;
    }

    for(k=0;k<ncil*j;k++)
    {
      aux  = (type_real)matrix_npart[i][k]*volinv[k/j];
      aux  -= 1.0f;

      matrix_mean_par[i][k] *= (matrix_npart[i][k]>0) ? 1./(type_real)matrix_npart[i][k] : 0.0f;
      matrix_mean_per[i][k] *= (matrix_npart[i][k]>0) ? 1./(type_real)matrix_npart[i][k] : 0.0f;

      matrix_quad_par[i][k] = matrix_npart[i][k]>10 ? sqrt((matrix_quad_par[i][k] - \
      matrix_npart[i][k]*matrix_mean_par[i][k]*matrix_mean_par[i][k])*(1.f/(type_real)(matrix_npart[i][k]-1))) : 0.0f;

      matrix_quad_per[i][k] = matrix_npart[i][k]>10 ? sqrt((matrix_quad_per[i][k] - \
      matrix_npart[i][k]*matrix_mean_per[i][k]*matrix_mean_per[i][k])*(1.f/(type_real)(matrix_npart[i][k]-1))) : 0.0f;

      fwrite(&matrix_npart[i][k],sizeof(type_int),1,pf);
      fwrite(&aux,sizeof(type_real),1,pf);
      fwrite(&matrix_mean_par[i][k],sizeof(type_real),1,pf);
      fwrite(&matrix_quad_par[i][k],sizeof(type_real),1,pf);
      fwrite(&matrix_mean_per[i][k],sizeof(type_real),1,pf);
      fwrite(&matrix_quad_per[i][k],sizeof(type_real),1,pf);
    }

    free(matrix_mean_par[i]);
    free(matrix_quad_par[i]);
    free(matrix_mean_per[i]);
    free(matrix_quad_per[i]);
    free(matrix_npart[i]);

    #ifdef LOCK
      omp_destroy_lock(&(lock[i]));
    #endif
  }

  fclose(pf);

  free(matrix_mean_per);
  free(matrix_quad_per);
  free(matrix_mean_par);
  free(matrix_quad_par);
  free(matrix_npart);
  free(matrix_middle);
  #ifdef LOCK
    free(lock);
  #endif

  free(rcil2);
  free(volinv);

  return;
}
