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

static type_int cpy;
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

static void set_name(char * name, const type_int NNN, const char * prefix)
{

  sprintf(name,"%.2d_%.4d",snap.num,NNN);

  sprintf(name,"%s_type_%.2d_pares",name,TYPE_FLAG);

  #ifdef MCRITIC
    sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%s",name,prefix);

  #ifdef CALCULA_VCM
    sprintf(name,"%s_vcm",name);
  #endif

  #ifdef EXTEND
    sprintf(name,"%s_extend",name);
  #endif

  #ifdef FIXED_SEPARATION
    sprintf(name,"%s_%.2f_%.2f",name,RLEN,RSEP);
  #else
    sprintf(name,"%s_%.2d_%.2d",name,nbins,ncil);
  #endif

  sprintf(name,"%s_%.2f_Mpc",name,RAUX/1000.f);

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

static void sum_matrix(const type_int idpar,const type_int idfil, const type_int start, const type_int end, const type_real * restrict rcil2, type_real vdir[], const type_real mod_r)
{
  int k, j, idim, ix, iy, iseg;
  type_real racum_min,rx_min,ry_min;
  type_real r,rsep,rlen,racum,dot,dis;
  type_real Pos_cent[3], vprima[3], delta[3];

  iseg = -1;
  dis = racum = 0.0f;
  ry_min = cp.lbox*cp.lbox;
  racum_min = rx_min = -cp.lbox;
  rlen = mod_r;

  #ifdef EXTEND
  rlen += (mod_r);
  #endif

  ////////////////////////////////////////////////////////////////////////////////////////////

  //
  //REGION EXTEND
  //
  
  //rsep   = RAUX;
  rsep = 0.5*mod_r;

  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    Pos_cent[idim] = Seg_cent[start].Pos[idim] - rsep*vdir[idim]; // cambio el origen
    #ifdef PERIODIC
    Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
    Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
    #endif
  
    delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
    #ifdef PERIODIC
    if(delta[idim]> 0.5*cp.lbox) delta[idim] -= cp.lbox;
    if(delta[idim]<-0.5*cp.lbox) delta[idim] += cp.lbox;
    #endif
    r += delta[idim]*delta[idim];
  }

  dot = delta[0]*vdir[0]+delta[1]*vdir[1]+delta[2]*vdir[2];

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

  #ifdef EXTEND
  racum += rsep;
  #endif

  ////////////////////////////////////////////////////////////////////////////////////////////

  //
  // REGION INTERNA
  //

  rsep = mod_r;
 
  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    Pos_cent[idim] = Seg_cent[start].Pos[idim];
    #ifdef PERIODIC
    Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
    Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
    #endif

    delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
    #ifdef PERIODIC
    if(delta[idim]> 0.5*cp.lbox) delta[idim] -= cp.lbox;
    if(delta[idim]<-0.5*cp.lbox) delta[idim] += cp.lbox;
    #endif
    r += delta[idim]*delta[idim];
  }

  dot = delta[0]*vdir[0]+delta[1]*vdir[1]+delta[2]*vdir[2];

  // SI ENTRA EN MI CILINDRO
  if(dot>=0.0f && dot<=rsep) 
  {
    dis = fabs(r - dot*dot);
    if(dis<ry_min && dis<rcil2[0])
    {
      ry_min = dis;
      rx_min = dot; 
      racum_min = racum; 
      iseg = 1;
    }
  }

  racum += rsep;
  

  ////////////////////////////////////////////////////////////////////////////////////////////
  
  //
  //REGION EXTEND
  //

  //rsep   = RAUX;
  rsep = 0.5*mod_r;

  r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    Pos_cent[idim] = Seg_cent[end].Pos[idim];
    #ifdef PERIODIC
    Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
    Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
    #endif

    delta[idim] = P[idpar].Pos[idim] - Pos_cent[idim];
    #ifdef PERIODIC
    if(delta[idim]> 0.5*cp.lbox) delta[idim] -= cp.lbox;
    if(delta[idim]<-0.5*cp.lbox) delta[idim] += cp.lbox;
    #endif
    r += delta[idim]*delta[idim];
  }

  dot = delta[0]*vdir[0]+delta[1]*vdir[1]+delta[2]*vdir[2];

  // SI ENTRA EN MI CILINDRO
  if(dot>=0.0f && dot<=rsep) 
  {
    dis = fabs(r - dot*dot);
    if(dis<ry_min && dis<rcil2[0])
    {
      ry_min = dis;
      rx_min = dot;
      racum_min = racum; 
      iseg = 2;
    }
  }

  #ifdef EXTEND
  racum += rsep;
  #endif

  if(iseg<0)
    return;

  j = nbins; 

  #ifdef EXTEND
    j *= 2;
  #endif

  if(iseg==0)
  {
    #ifndef EXTEND
      return;
    #endif
   
    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] = Seg_cent[start].Pos[idim] - RAUX*vdir[idim]; // cambio el origen

      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

  }else if(iseg==2){

    #ifndef EXTEND
      return;
    #endif

    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] = Seg_cent[end].Pos[idim]; 

      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

  }else{

    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] = Seg_cent[start].Pos[idim]; 

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
   //if(cpy!=0)
   //  vdir[idim] *= -1.0f;

   vprima[idim] = P[idpar].Vel[idim]; // asigna

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
  if((iseg == 0) || (iseg == 2))
  {
    iseg = -1;
  }else
  #endif 
  {

    #ifdef EXTEND
      r = rsep-mod_r;
    #else
      r = rsep;
    #endif
    r *= (1.0/mod_r);

    if(r>0.25 && r<0.75)
    {
      if(r_mean<sqrt(ry_min))
        iseg = -1;
      else
        iseg =  1;
    }else{
      iseg = -1;
    }
    
  }

  k = idfil+cpy*cp.nseg;

  #ifdef LOCK
  omp_set_lock(&lock[idfil]);
  #endif
  matrix_mean_par[k][j] += dot;
  matrix_quad_par[k][j] += (dot*dot);
  matrix_mean_per[k][j] += racum;
  matrix_quad_per[k][j] += (racum*racum);  
  matrix_npart[k][j]++;
  if(iseg>0) 
    matrix_middle[k]++;
  #ifdef LOCK
  omp_unset_lock(&lock[idfil]);
  #endif

  return;
}

static void search_vdir(type_real vdir[], type_real *r, type_int init, type_int end)
{
  type_int  idim;

  *r = 0.0f;
  for(idim=0;idim<3;idim++)
  {
    vdir[idim] = Seg_cent[end].Pos[idim]-Seg_cent[init].Pos[idim];    
    #ifdef PERIODIC
    if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
    if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
    #endif
    *r += vdir[idim]*vdir[idim];
  }
  *r = sqrt(*r);

  for(idim=0;idim<3;idim++)
    vdir[idim] /= (*r); 

  return;

}

static void calc_fil_dis(const type_int idpar, const type_int idnod, const type_real * restrict rcil2)
{
  const type_int idfil = idnod;
  type_int n,start,end;
  type_real vdir[3], r;

  if(cpy==0)
  {

    for(n=1;n<=NRAND;n++)
    {
      start = idnod;
      end   = idnod+n*cp.nseg;
      search_vdir(vdir,&r,start,end);    
      sum_matrix(idpar,idfil,start,end,rcil2,vdir,r);
    }

  }else{

    for(n=1;n<=NRAND;n++)
    {
      start = idnod+n*cp.nseg;
      end   = idnod;
      search_vdir(vdir,&r,start,end);
      sum_matrix(idpar,idfil,start,end,rcil2,vdir,r);
    }

  }
  
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
          Posprima[0] = Seg_cent[i].Pos[0] - Pos_cent[0];
          Posprima[1] = Seg_cent[i].Pos[1] - Pos_cent[1];
          Posprima[2] = Seg_cent[i].Pos[2] - Pos_cent[2];

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
            //dat.id = (i%cp.nseg);
            dat.id = Seg_cent[i].Id;
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

static void calcular_pert(const type_int i, const type_real rbus)
{
  struct list *root;
  struct list *aux_lista_fil = NULL;
  type_real *rcil2;
  type_real r;

  calc_search_fil(P[i].Pos,rbus*rbus,&aux_lista_fil);             
  remove_duplicates(aux_lista_fil);

  rcil2  = (type_real *) malloc((ncil+1)*sizeof(type_real));

  r = 0.5*rbus;

  #ifdef BIN_LOG
    logspace(rcil2,r,r/50.0f,ncil+1);
  #else
    linspace(rcil2,r,0.0f,ncil+1);
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

static void build_rand(void)
{
  type_int i,j,k,l,idim;
  type_real alfa, delta;
  type_real vdir[3], r;
  type_int  count = 0;

  for(i=0;i<cp.nseg;i++)
  { 

    Seg_cent[i].Id = i;
    memcpy(Seg_cent[i].Pos,Gr[Seg[i].list[cpy]].Pos,3*sizeof(type_real)); 

    for(j=1;j<=NRAND;j++)
    {
      alfa  = drand48()*2.0f*M_PI;
  	  delta = 2.0*drand48()-1.0f;
	    delta = asin(delta);

      vdir[0] = cos(delta)*cos(alfa);
      vdir[1] = cos(delta)*sin(alfa);
      vdir[2] = sin(delta);    

      k = i+j*cp.nseg;
      Seg_cent[k].Id = i;

      for(idim=0;idim<3;idim++)
      {
        r = Seg[i].len;
        Seg_cent[k].Pos[idim] = Seg_cent[i].Pos[idim] + r*vdir[idim];
        #ifdef PERIODIC
        Seg_cent[k].Pos[idim] = Seg_cent[k].Pos[idim] >= cp.lbox ? Seg_cent[k].Pos[idim]-cp.lbox : Seg_cent[k].Pos[idim];
        Seg_cent[k].Pos[idim] = Seg_cent[k].Pos[idim] <      0.0 ? Seg_cent[k].Pos[idim]+cp.lbox : Seg_cent[k].Pos[idim];
        #endif
      }
    }
    
    if(Seg[i].len>2.0*RAUX)
    { 
      j = (type_int)((0.5*Seg[i].len)/RAUX); 
      j = j == 0 ? 2*NRAND : 2*j*NRAND;
      count += j;
    }
  }

  cp.ncen_rand += count;
  Seg_cent = (struct centr_data *) realloc(Seg_cent,cp.ncen_rand*sizeof(struct centr_data));
  
  count = (NRAND+1)*cp.nseg;

  for(i=0;i<cp.nseg;i++)
  { 
    if(Seg[i].len>2.0*RAUX)
    { 
      k  = (type_int)((0.5*Seg[i].len)/RAUX);
      k = k == 0 ? 1 : k;
      
      if(cpy==0)
      {

        for(j=1;j<=NRAND;j++)
        {
          search_vdir(vdir,&r,i,i+j*cp.nseg);    

          for(l=0;l<k;l++)
          {
            r = (type_real)(l+1)*RAUX;
            for(idim=0;idim<3;idim++)
            {
              Seg_cent[count].Pos[idim] = Seg_cent[i].Pos[idim] - r*vdir[idim];
              #ifdef PERIODIC
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] >= cp.lbox ? Seg_cent[count].Pos[idim]-cp.lbox : Seg_cent[count].Pos[idim];
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] <      0.0 ? Seg_cent[count].Pos[idim]+cp.lbox : Seg_cent[count].Pos[idim];
              #endif
            }
            Seg_cent[count].Id = i;
            count++;

            for(idim=0;idim<3;idim++)
            {
              Seg_cent[count].Pos[idim] = Seg_cent[i+j*cp.nseg].Pos[idim] + r*vdir[idim];
              #ifdef PERIODIC
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] >= cp.lbox ? Seg_cent[count].Pos[idim]-cp.lbox : Seg_cent[count].Pos[idim];
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] <      0.0 ? Seg_cent[count].Pos[idim]+cp.lbox : Seg_cent[count].Pos[idim];
              #endif
            }
            Seg_cent[count].Id = i;
            count++;
          }
        }

      }else{

        for(j=1;j<=NRAND;j++)
        {
          search_vdir(vdir,&r,i+j*cp.nseg,i);

          for(l=0;l<k;l++)
          {
            r = (type_real)(l+1)*RAUX;
            for(idim=0;idim<3;idim++)
            {
              Seg_cent[count].Pos[idim] = Seg_cent[i].Pos[idim] - r*vdir[idim];
              #ifdef PERIODIC
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] >= cp.lbox ? Seg_cent[count].Pos[idim]-cp.lbox : Seg_cent[count].Pos[idim];
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] <      0.0 ? Seg_cent[count].Pos[idim]+cp.lbox : Seg_cent[count].Pos[idim];
              #endif
            }
            Seg_cent[count].Id = i;
            count++;

            for(idim=0;idim<3;idim++)
            {
              Seg_cent[count].Pos[idim] = Seg_cent[i+j*cp.nseg].Pos[idim] + r*vdir[idim];
              #ifdef PERIODIC
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] >= cp.lbox ? Seg_cent[count].Pos[idim]-cp.lbox : Seg_cent[count].Pos[idim];
              Seg_cent[count].Pos[idim] = Seg_cent[count].Pos[idim] <      0.0 ? Seg_cent[count].Pos[idim]+cp.lbox : Seg_cent[count].Pos[idim];
              #endif
            }
            Seg_cent[count].Id = i;
            count++;
          }
        }
      }
    }

  }






  return;

}

static void init_rand(const type_int orden, const type_real aux)
{
  type_int i;
  #ifdef CALCULA_VCM
  type_int k;
  #endif

  cpy = orden;

  GREEN("********** IMPORTANTE ***********\n");
  sprintf(message,"***** CREA RANDOM FIJO EL %d *****\n",cpy);
  RED(message);
  
  build_rand();

  GREEN("**********************************\n");

  grid.nobj = cp.ncen_rand;
  grid.ngrid = (long)(cp.lbox/aux);

  if(grid.ngrid > NGRIDMAX)
  {
    grid.ngrid = NGRIDMAX;
    fprintf(stdout,"Using NGRIDMAX = %lu\n",grid.ngrid);
  }else{
    fprintf(stdout,"Using NGRID    = %lu\n",grid.ngrid);
  }
  fflush(stdout);

  grid_init();
  grid_build();

  BLUE("******************************\n");
  GREEN("******************************\n");
  fflush(stdout);

  #ifdef CALCULA_VCM

    GREEN("CALCULA VCM\n");

    for(i=0;i<cp.nseg;i++)  
    {
      for(k=0;k<3;k++)
      {
        //Seg[i].Vmedia[k] = (Seg[i].Mass[0]*Seg[i].Vnodos[k]+    \
        //                    Seg[i].Mass[1]*Seg[i].Vnodos[k+3])/ \
        //                   (Seg[i].Mass[0]+Seg[i].Mass[1]);
        
        Seg[i].Vmedia[k] = Seg[i].Vnodos[k+3*cpy];

      }
    }

    GREEN("Termina el calculo\n");

  #endif

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) \
  private(i) shared(cp,stdout)
  for(i=0;i<cp.npart;i++)
  {
    if(i%1000000==0) fprintf(stdout,"%u %.4f\n",i,(float)i/(float)cp.npart);
    calcular_pert(i,aux);
  }

  grid_free();

  GREEN("*********************************\n");
  sprintf(message,"****** TERMINO EL CALCULO %d *****\n",cpy);
  RED(message);
  GREEN("**********************************\n");
  fflush(stdout);

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
  #ifdef NODOS
    FILE *pf_left, *pf_right;
  #endif

  srand48(80); // guarda la semilla

  j = nbins;

  #ifdef EXTEND
    j *= 2;  
  #endif

  cp.ncen_rand = (NRAND+1)*cp.nseg;

  Seg_cent = (struct centr_data *) malloc(cp.ncen_rand*sizeof(struct centr_data));
  matrix_mean_per = (type_real **) malloc(2*cp.nseg*sizeof(type_real *));
  matrix_quad_per = (type_real **) malloc(2*cp.nseg*sizeof(type_real *));
  matrix_mean_par = (type_real **) malloc(2*cp.nseg*sizeof(type_real *));
  matrix_quad_par = (type_real **) malloc(2*cp.nseg*sizeof(type_real *));
  matrix_npart    = (type_int  **) malloc(2*cp.nseg*sizeof(type_int  *));
  matrix_middle   = (type_int  *)  calloc(2*cp.nseg,sizeof(type_int));
  #ifdef LOCK
    lock = (omp_lock_t *) malloc(cp.nseg*sizeof(omp_lock_t));
  #endif

  for(i=0;i<cp.nseg;i++)  
  {
    matrix_mean_per[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_mean_per[i+cp.nseg] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_quad_per[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_quad_per[i+cp.nseg] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_mean_par[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_mean_par[i+cp.nseg] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_quad_par[i] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_quad_par[i+cp.nseg] = (type_real *) calloc(ncil*j,sizeof(type_real));
    matrix_npart[i]    = (type_int  *) calloc(ncil*j,sizeof(type_int ));
    matrix_npart[i+cp.nseg]    = (type_int  *) calloc(ncil*j,sizeof(type_int ));
    #ifdef LOCK
      omp_init_lock(&(lock[i]));
    #endif
  }

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  BLUE("*********************\n");
  init_rand(0,r_aux);
  BLUE("*********************\n");
  init_rand(1,r_aux);
  BLUE("*********************\n");

  GREEN("Comienza a escribir\n");

  rcil2  = (type_real *) malloc((ncil+1)*sizeof(type_real));
  volinv = (type_real *) malloc(ncil*sizeof(type_real));

  r = RAUX; // RADIO DE BUSQUEDA EN Kpc
  #ifdef BIN_LOG
    logspace(rcil2,r,r/50.0f,ncil+1);
  #else
    linspace(rcil2,r,0.0f,ncil+1);
  #endif


  #ifdef NODOS

    set_name(filename,NNN,"matrix_left");
    pf_left = fopen(filename,"w");

    fwrite(&cp.nseg,sizeof(type_int),1,pf_left);        
    fwrite(&j,sizeof(type_int),1,pf_left);        
    fwrite(&ncil,sizeof(type_int),1,pf_left);        

    set_name(filename,NNN,"matrix_right");
    pf_right = fopen(filename,"w");

    fwrite(&cp.nseg,sizeof(type_int),1,pf_right);
    fwrite(&j,sizeof(type_int),1,pf_right);      
    fwrite(&ncil,sizeof(type_int),1,pf_right);   

  #endif

  set_name(filename,NNN,"matrix");
  pf = fopen(filename,"w");

  fwrite(&cp.nseg,sizeof(type_int),1,pf);        
  fwrite(&j,sizeof(type_int),1,pf);        
  fwrite(&ncil,sizeof(type_int),1,pf);        

  for(i=0;i<cp.nseg;i++) 
  {
    #ifdef NODOS

      fwrite(&i,sizeof(type_int),1,pf_left);            // Num Filamento left
      fwrite(&Seg[i].flag,sizeof(type_int),1,pf_left);  // escribe la bandera left
      fwrite(&Seg[i].len,sizeof(type_real),1,pf_left);  // escribe la longitud del segmento left

      fwrite(&i,sizeof(type_int),1,pf_right);           // Num Filamento right 
      fwrite(&Seg[i].flag,sizeof(type_int),1,pf_right); // escribe la bandera right 
      fwrite(&Seg[i].len,sizeof(type_real),1,pf_right); // escribe la longitud del segmento right 

    #endif

    fwrite(&i,sizeof(type_int),1,pf);              // Num Filamento
    fwrite(&Seg[i].flag,sizeof(type_int),1,pf);    // escribe la bandera
    fwrite(&Seg[i].len,sizeof(type_real),1,pf);    // escribe la longitud del segmento

    #ifdef NODOS

      fwrite(&Seg[i].Mass[0],sizeof(type_real),1,pf_left);    // escribe la masa del primer nodo left 
      fwrite(&Seg[i].Mass[1],sizeof(type_real),1,pf_left);    // escribe la masa del ultimo nodo left   

      fwrite(&Seg[i].Mass[0],sizeof(type_real),1,pf_right);   // escribe la masa del primer nodo right 
      fwrite(&Seg[i].Mass[1],sizeof(type_real),1,pf_right);   // escribe la masa del ultimo nodo right   

    #endif
   
    fwrite(&Seg[i].Mass[0],sizeof(type_real),1,pf);         // escribe la masa del primer nodo
    fwrite(&Seg[i].Mass[1],sizeof(type_real),1,pf);         // escribe la masa del ultimo nodo   

    matrix_middle[i]         /= NRAND;
    matrix_middle[i+cp.nseg] /= NRAND;
    
    #ifdef NODOS

      aux  = (type_real)matrix_middle[i]*(cp.lbox*cp.lbox*cp.lbox);
      aux /= (type_real)cp.npart;
      aux /= (Seg[i].len*0.5*4.0*M_PI*r_mean*r_mean);
      aux -= 1.0f;
      fwrite(&aux,sizeof(type_real),1,pf_left);  // escribe la masa dentro de un cilindro de r_mean left

      aux  = (type_real)matrix_middle[i+cp.nseg]*(cp.lbox*cp.lbox*cp.lbox);
      aux /= (type_real)cp.npart;
      aux /= (Seg[i].len*0.5*4.0*M_PI*r_mean*r_mean);
      aux -= 1.0f;
      fwrite(&aux,sizeof(type_real),1,pf_right); // escribe la masa dentro de un cilindro de r_mean right

    #endif

    matrix_middle[i]         += matrix_middle[i+cp.nseg];

    aux  = (type_real)matrix_middle[i]*(cp.lbox*cp.lbox*cp.lbox);
    aux /= (type_real)cp.npart;
    aux /= (Seg[i].len*0.5*4.0*M_PI*r_mean*r_mean);
    aux -= 1.0f;
    fwrite(&aux,sizeof(type_real),1,pf);             // escribe la masa dentro de un cilindro de r_mean

    rlong = Seg[i].len;
  
    #ifdef EXTEND
    rlong += (2.0f*RAUX);
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
      ///////////////////////////////////////////////

      matrix_npart[i][k]    /= NRAND;
      matrix_mean_par[i][k] /= NRAND;
      matrix_quad_par[i][k] /= NRAND;
      matrix_mean_per[i][k] /= NRAND;
      matrix_quad_per[i][k] /= NRAND;

      ///////////////////////////////////////////////
      ///////////////////////////////////////////////
     
      matrix_npart[i+cp.nseg][k]    /= NRAND;
      matrix_mean_par[i+cp.nseg][k] /= NRAND;
      matrix_quad_par[i+cp.nseg][k] /= NRAND;
      matrix_mean_per[i+cp.nseg][k] /= NRAND;
      matrix_quad_per[i+cp.nseg][k] /= NRAND;

      ///////////////////////////////////////////////
 
      #ifdef NODOS

        ///////////////////////////////////////////////////////////////////////////////////////

        fwrite(&matrix_npart[i][k],sizeof(type_int),1,pf_left); // npart

        aux  = (type_real)matrix_npart[i][k]*volinv[k/j];
        aux  -= 1.0f;

        fwrite(&aux,sizeof(type_real),1,pf_left); // rho

        aux  = matrix_mean_par[i][k];
        aux *= matrix_npart[i][k]>0 ? 1./(type_real)matrix_npart[i][k] : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_left); // matrix_mean_par

        aux = matrix_npart[i][k]>10 ? sqrt(fabs(matrix_quad_par[i][k] - \
        matrix_npart[i][k]*aux*aux)*(1.f/(type_real)(matrix_npart[i][k]-1))) : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_left); // matrix_quad_par

        aux  = matrix_mean_per[i][k];
        aux *= matrix_npart[i][k]>0 ? 1./(type_real)matrix_npart[i][k] : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_left); // matrix_mean_per

        aux = matrix_npart[i][k]>10 ? sqrt(fabs(matrix_quad_per[i][k] - \
        matrix_npart[i][k]*aux*aux)*(1.f/(type_real)(matrix_npart[i][k]-1))) : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_left); // matrix_quad_per 

        ///////////////////////////////////////////////////////////////////////////////////////

        fwrite(&matrix_npart[i+cp.nseg][k],sizeof(type_int),1,pf_right); // npart

        aux  = (type_real)matrix_npart[i+cp.nseg][k]*volinv[k/j];
        aux  -= 1.0f;

        fwrite(&aux,sizeof(type_real),1,pf_right); // rho

        aux  = matrix_mean_par[i+cp.nseg][k];
        aux *= matrix_npart[i+cp.nseg][k]>0 ? 1./(type_real)matrix_npart[i+cp.nseg][k] : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_right); // matrix_mean_par

        aux = matrix_npart[i+cp.nseg][k]>10 ? sqrt(fabs(matrix_quad_par[i+cp.nseg][k] - \
        matrix_npart[i+cp.nseg][k]*aux*aux)*(1.f/(type_real)(matrix_npart[i+cp.nseg][k]-1))) : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_right); // matrix_quad_par

        aux  = matrix_mean_per[i+cp.nseg][k];
        aux *= matrix_npart[i+cp.nseg][k]>0 ? 1./(type_real)matrix_npart[i+cp.nseg][k] : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_right); // matrix_mean_per

        aux = matrix_npart[i+cp.nseg][k]>10 ? sqrt(fabs(matrix_quad_per[i+cp.nseg][k] - \
        matrix_npart[i+cp.nseg][k]*aux*aux)*(1.f/(type_real)(matrix_npart[i+cp.nseg][k]-1))) : 0.0f;

        fwrite(&aux,sizeof(type_real),1,pf_right); // matrix_quad_per

        ///////////////////////////////////////////////////////////////////////////////////////
 
      #endif

      matrix_npart[i][k]    += matrix_npart[i+cp.nseg][k];
      matrix_mean_par[i][k] += matrix_mean_par[i+cp.nseg][k];
      matrix_quad_par[i][k] += matrix_quad_par[i+cp.nseg][k];
      matrix_mean_per[i][k] += matrix_mean_per[i+cp.nseg][k];
      matrix_quad_per[i][k] += matrix_quad_per[i+cp.nseg][k];
      
      //
      //

      matrix_mean_par[i][k] *= (matrix_npart[i][k]>0) ? 1./(type_real)matrix_npart[i][k] : 0.0f;
      matrix_mean_per[i][k] *= (matrix_npart[i][k]>0) ? 1./(type_real)matrix_npart[i][k] : 0.0f;

      matrix_quad_par[i][k] = matrix_npart[i][k]>10 ? sqrt(fabs(matrix_quad_par[i][k] - \
      matrix_npart[i][k]*matrix_mean_par[i][k]*matrix_mean_par[i][k])*(1.f/(type_real)(matrix_npart[i][k]-1))) : 0.0f;

      matrix_quad_per[i][k] = matrix_npart[i][k]>10 ? sqrt(fabs(matrix_quad_per[i][k] - \
      matrix_npart[i][k]*matrix_mean_per[i][k]*matrix_mean_per[i][k])*(1.f/(type_real)(matrix_npart[i][k]-1))) : 0.0f;

      //
      //

      matrix_npart[i][k]    *= 0.5; 
      matrix_mean_par[i][k] *= 0.5;
      matrix_quad_par[i][k] *= 0.5;
      matrix_mean_per[i][k] *= 0.5;
      matrix_quad_per[i][k] *= 0.5;

      aux  = (type_real)matrix_npart[i][k]*volinv[k/j];
      aux  -= 1.0f;

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

    free(matrix_mean_par[i+cp.nseg]);
    free(matrix_quad_par[i+cp.nseg]);
    free(matrix_mean_per[i+cp.nseg]);
    free(matrix_quad_per[i+cp.nseg]);
    free(matrix_npart[i+cp.nseg]);

    #ifdef LOCK
      omp_destroy_lock(&(lock[i]));
    #endif
  }

  fclose(pf);
  #ifdef NODOS
    fclose(pf_left);
    fclose(pf_right);
  #endif

  GREEN("Termina de escribir\n");
  BLUE("*********************\n");

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
