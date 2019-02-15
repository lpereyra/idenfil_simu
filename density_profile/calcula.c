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
#include "timer.h"

#include <sys/stat.h>
#include <sys/types.h>

#define IX(i,j,n) n*j+i

//////////////////////////////
static type_int  **matrix_npart;
//////////////////////////////
#ifdef LOCK
  static omp_lock_t *lock; 
#endif

static void set_name(char * name, const type_int NNN, const char * prefix)
{
  struct stat st = {0};

  #ifdef RANDOM
    sprintf(name,"./random/");
  #else
    sprintf(name,"./data/");
  #endif

  if(stat(name, &st) == -1) 
  {
    mkdir(name, 0700);
  }

  #ifdef ORIGINAL
    sprintf(name,"%s%.2d_%.4d_type_%.2d",name,snap.num,NNN,TYPE_FLAG);
  #else
    sprintf(name,"%s%.2d_%.4d_type_%.2d_smooth",name,snap.num,NNN,TYPE_FLAG);
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

  #ifdef EXTEND
    sprintf(name,"%s_extend",name);
  #endif

  #ifdef RANDOM
    sprintf(name,"%s_%.2d_random",name,NRANDOM);
  #else
    sprintf(name,"%s_data",name);
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

static void calc_fil_dis(const type_real * restrict Pos_idpar, const type_int idfil, const type_real * restrict rcil2)
{
  int j, k, idim, ix, iy, iseg;
  type_real racum_min,rx_min,ry_min;
  type_real r,rsep,rlen,racum,dot,dis;
  type_real Pos_cent[3], vdir[3], delta[3];

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
    delta[idim] = Pos_idpar[idim] - Pos_cent[idim];
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
      delta[idim] = Pos_idpar[idim] - Pos_cent[idim];
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
    delta[idim] = Pos_idpar[idim] - Pos_cent[idim];
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

  #ifdef LOCK
  omp_set_lock(&lock[idfil]);
  #endif
  matrix_npart[idfil][j]++;
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

static void calcular_pert(const type_real * restrict Pos_idpar, const type_real rbus)
{
  struct list *root;
  struct list *aux_lista_fil = NULL;
  type_real *rcil2;

  calc_search_fil(Pos_idpar,rbus*rbus,&aux_lista_fil);             

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
    calc_fil_dis(Pos_idpar,root->data.id,rcil2);
    aux_lista_fil = aux_lista_fil->next;
    free(root);
  }

  free(rcil2);

  return;
}

extern void propiedades(type_int NNN, type_real *fof)
{
  long i,j,k;
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

  matrix_npart    = (type_int  **) malloc(cp.nseg*sizeof(type_int  *));
  #ifdef LOCK
    lock = (omp_lock_t *) malloc(cp.nseg*sizeof(omp_lock_t));
  #endif
  aux  = cbrt(3.0/(4.0*M_PI*cp.Mpart*(type_real)cp.npart));
  aux *= cp.lbox;
  aux *= fof[1];

  for(i=0;i<cp.nseg;i++)  
  {
    matrix_npart[i]    = (type_int  *) calloc(ncil*j,sizeof(type_int ));
    Seg[i].Rvir[0] = cbrt(Seg[i].Mass[0])*aux;
    Seg[i].Rvir[1] = cbrt(Seg[i].Mass[1])*aux;
    Seg[i].Rvir[0] *= Seg[i].Rvir[0];
    Seg[i].Rvir[1] *= Seg[i].Rvir[1];
    #ifdef LOCK
      omp_init_lock(&(lock[i]));
    #endif
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

  #ifdef RANDOM

    long *nrand;
    nrand = (long *) malloc(NTHREADS*sizeof(long));

    pf=fopen("/dev/urandom","r");
    for(i=0;i<NTHREADS;i++)
      fread(&nrand[i],sizeof(long),1,pf);
    fclose(pf);

    #ifdef EXTEND
      sprintf(filename,"seed_extend_type_%.2d.dat",TYPE_FLAG);
    #else
      sprintf(filename,"seed_type_%.2d.dat",TYPE_FLAG);
    #endif
    pf = fopen(filename,"w");

    fprintf(pf,"%d\n",NTHREADS);
    for(i=0;i<NTHREADS;i++)
      fprintf(pf,"%d %ld\n",(int)i,nrand[i]);
    fclose(pf);
  
    k = 0;
    #pragma omp parallel \
    num_threads(NTHREADS) default(none) \
    private(i) shared(cp,nrand,stdout) \
    reduction(+:k)
    {
      type_int  Tid = omp_get_thread_num(); 
      type_real Pos_rand[3];

      srand48(nrand[Tid]); // guarda la semilla
    
      nrand[Tid] = NRANDOM*(long)cp.npart/NTHREADS;
      nrand[Tid] = Tid==NTHREADS-1 ? NRANDOM*(long)cp.npart-(NTHREADS-1)*nrand[Tid] : nrand[Tid];

      for(i=0;i<nrand[Tid];i++)
      {
        if(i%10000000==0) fprintf(stdout,"%u %.4f\n",Tid,(float)i/(float)nrand[Tid]);
        Pos_rand[0] = drand48()*cp.lbox;
        Pos_rand[1] = drand48()*cp.lbox;
        Pos_rand[2] = drand48()*cp.lbox;
        calcular_pert(Pos_rand,r_aux);
        k++;
      }

    }
    
    assert((NRANDOM*(long)cp.npart)==k);
    free(nrand);

  #else

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) \
    private(i) shared(cp,P,stdout)
    for(i=0;i<cp.npart;i++)
    {
      if(i%1000000==0) fprintf(stdout,"%ld %.4f\n",i,(float)i/(float)cp.npart);
      calcular_pert(P[i].Pos,r_aux);
    }

  #endif

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
    Seg[i].Rvir[0] = sqrt(Seg[i].Rvir[0]);
    Seg[i].Rvir[1] = sqrt(Seg[i].Rvir[1]);
    
    ////////////////////////////////////////////////////////////////

    fwrite(&i,sizeof(type_int),1,pf);                   // Num Filamento
    fwrite(&Seg[i].flag,sizeof(type_int),1,pf);         // escribe la bandera
    fwrite(&Seg[i].len,sizeof(type_real),1,pf);         // escribe la longitud del segmento

    fwrite(&Seg[i].Mass[0],sizeof(type_real),1,pf);             // escribe la masa del primer nodo
    fwrite(&Seg[i].Mass[1],sizeof(type_real),1,pf);             // escribe la masa del ultimo nodo   

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
        //if(k != ncil-1) // para hacer el ultimo cilindro entero
        volinv[k] -= r;
      #endif  

      volinv[k] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[k])/(type_real)cp.npart;
    }

    for(k=0;k<ncil*j;k++)
    {
      aux  = (type_real)matrix_npart[i][k]*volinv[k/j];
      aux  -= 1.0f;

      fwrite(&matrix_npart[i][k],sizeof(type_int),1,pf);
      fwrite(&aux,sizeof(type_real),1,pf);
    }

    free(matrix_npart[i]);

    #ifdef LOCK
      omp_destroy_lock(&(lock[i]));
    #endif
  }

  fclose(pf);

  free(matrix_npart);
  #ifdef LOCK
    free(lock);
  #endif

  free(rcil2);
  free(volinv);

  return;
}
