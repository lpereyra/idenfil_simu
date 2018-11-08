#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "cosmoparam.h"
#include "variables.h"
#include "leesnap.h"
#include "grid.h"
#include "iden.h"
#include "bitmask.h"

#ifdef LOCK
  static omp_lock_t *lock;
#endif

static inline long igrid(const long ix, const long iy, const long iz, const long ngrid) 
{
  return (ix*ngrid+ iy)*ngrid + iz;
}

static inline type_int Raiz(type_int i, type_int * restrict ar)
{

 if(i != ar[i])
   ar[i] = Raiz(ar[i],ar);

 return ar[i];

}

static inline void Unir(type_int u, type_int v, type_int * restrict ar)
{
 
  type_int z;

  while(ar[u] != ar[v])
  { 
      if(ar[u] < ar[v])
      {
#ifdef LOCK
          if(u == ar[u])
          {
            omp_set_lock(&(lock[u]));
            z = 0;

            if(u == ar[u])
            {
              ar[u] = ar[v];  
              z = 1;
            } 

            omp_unset_lock(&(lock[u]));
            if(z==1) break;             

          }
#else
          if(u == ar[u])
          {
            ar[u] = ar[v];  
            break;             
          }
#endif
          
          z = ar[u];   
          ar[u] = ar[v];
          u = z;

      }else{
#ifdef LOCK
          if(v == ar[v])
          {
            omp_set_lock(&(lock[v]));
            z = 0;

            if(v == ar[v])
            {
              ar[v] = ar[u];   
              z = 1;
            }

            omp_unset_lock(&(lock[v]));
            if(z == 1) break;            
          }
#else
          if(v == ar[v])
          {
            ar[v] = ar[u];   
            break;             
          }
#endif

          z = ar[v];   
          ar[v] = ar[u];   
          v = z;

      }
  }

}

static void busv(const type_int ic, type_int * restrict test)
{

  long ixc, iyc, izc, ibox;
  type_int i, niv;
  type_real xx, yy, zz;

  ixc  = (long)(P.x[ic]*(type_real)grid.ngrid*(1.f/cp.lbox));
  iyc  = (long)(P.y[ic]*(type_real)grid.ngrid*(1.f/cp.lbox));
  izc  = (long)(P.z[ic]*(type_real)grid.ngrid*(1.f/cp.lbox));

  #ifndef PERIODIC
    for(long ixx = ((ixc-1<0) ? 0 : ixc-1); ixx <= ((ixc+1 >= grid.ngrid) ? grid.ngrid-1 : ixc+1); ixx++)
  #else
    for(long ixx = ixc-1; ixx <= ixc+1; ixx++)
  #endif
  {
    #ifndef PERIODIC
      for(long iyy = ((iyc-1<0) ? 0 : iyc-1); iyy <= ((iyc+1 >= grid.ngrid) ? grid.ngrid-1 : iyc+1); iyy++)
    #else
      for(long iyy = iyc-1 ; iyy <= iyc+1 ; iyy++)
    #endif
    {
      #ifndef PERIODIC
        for(long izz = ((izc-1<0) ? 0 : izc-1); izz <= ((izc+1 >= grid.ngrid) ? grid.ngrid-1 : izc+1); izz++)
      #else
        for(long izz = izc-1 ; izz <= izc+1 ; izz++)
      #endif
      {

      	#ifdef PERIODIC
          ibox = igrid(( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) ),\
                       ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ),\
                       ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) ),\
                       (long)grid.ngrid);
        #else
          ibox = igrid(ixx,iyy,izz,(long)grid.ngrid);
        #endif

        for(i=grid.icell[ibox];i<grid.icell[ibox]+grid.size[ibox];i++)
        {
          if(!TestBit(test,i))
          {

            xx = P.x[i] - P.x[ic];
            yy = P.y[i] - P.y[ic];
            zz = P.z[i] - P.z[ic];

            #ifdef PERIODIC
            xx = ( xx >  cp.lbox*0.5f ) ? xx - cp.lbox : xx ;
            yy = ( yy >  cp.lbox*0.5f ) ? yy - cp.lbox : yy ;
            zz = ( zz >  cp.lbox*0.5f ) ? zz - cp.lbox : zz ;
            xx = ( xx < -cp.lbox*0.5f ) ? xx + cp.lbox : xx ;
            yy = ( yy < -cp.lbox*0.5f ) ? yy + cp.lbox : yy ;
            zz = ( zz < -cp.lbox*0.5f ) ? zz + cp.lbox : zz ;
            #endif

            for(niv=0;niv<nfrac;niv++)
            {
      	      if(xx*xx + yy*yy + zz*zz < iden.r0[niv])
              {
                Unir(ic,i,gr[niv]);
              }
            }
            
          } // cierra el if

        } /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/

}

static void linkedlist(type_int *test, type_int * restrict ar)
{
  type_int i, g;

  Temp.ll = (type_int *) calloc(cp.npart,sizeof(type_int));

  iden.ngrupos = 0;
  for(i=0;i<cp.npart;i++)
  {
    ar[i] = Raiz(i,ar);
    if(TestBit(test,ar[i]))
    {
      Temp.ll[ar[i]]++;
      if(Temp.ll[ar[i]]>=NPARTMIN)
      { 
        iden.ngrupos++;
        Temp.ll[ar[i]] = iden.ngrupos;
        ClearBit(test,ar[i]);
      }
    }
  }

  iden.ngrupos++;  // SUMA UNO;

  Temp.head   = (type_int *) malloc(iden.ngrupos*sizeof(type_int));
  Temp.npgrup = (type_int *) malloc(iden.ngrupos*sizeof(type_int));

  for(i=0;i<iden.ngrupos;i++)
  {
    Temp.head[i]   = cp.npart;
    Temp.npgrup[i] =  0;
  }

  for(i=0;i<cp.npart;i++)
  { 
    if(TestBit(test,ar[i]))
    {
      ar[i] = 0;
    }else{
      ar[i] = Temp.ll[ar[i]];
    }
  }

  for(i=0;i<cp.npart;i++)
  {
    g = ar[i];

    SetBit(test,i);

    #ifdef DEBUG
    assert((g >= 0) && (g < iden.ngrupos));
    #endif
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
    Temp.npgrup[g]++;
  }

}

static void Write_Groups(type_int niv)
{
  type_int i,j,k,npar,gn,flag;
  type_real boxsize;
  char filename[200];
  FILE *pfout;

  ///////////////////////////////////////////////////////
  sprintf(filename,"%.2d_%.2f_fof.bin",snap.num,fof[niv]);
  pfout=fopen(filename,"w");
  //////////////////////////////////////////////////////

  j = 0;
  for(i=1;i<iden.ngrupos;i++)
    j +=  Temp.npgrup[i];

  header.npart[1] = j;
  header.npartTotal[1] = j;
  header.num_files = 1;

  boxsize = pmax[0] - pmin[0];
  if((pmax[1] - pmin[1]) > boxsize)
    boxsize = pmax[1] - pmin[1];
  if((pmax[2] - pmin[2]) > boxsize)
    boxsize = pmax[2] - pmin[2];

  header.BoxSize = boxsize;

  flag = sizeof(header);
  fwrite(&flag,sizeof(flag),1,pfout);
  fwrite(&header,sizeof(header),1,pfout);
  fwrite(&flag,sizeof(flag),1,pfout);

  flag = 3*header.npart[1]*sizeof(type_real);
  fwrite(&flag,sizeof(flag),1,pfout);

  npar = gn = 0;
  for(i=1;i<iden.ngrupos;i++)
  {
    j = 0;
    k = Temp.head[i];

    while(k != cp.npart)
    {
      fwrite(&P.x[k],sizeof(type_real),1,pfout);
      fwrite(&P.y[k],sizeof(type_real),1,pfout);
      fwrite(&P.z[k],sizeof(type_real),1,pfout);
      k = Temp.ll[k];
      j++;
    }

    assert(j == Temp.npgrup[i]);

    npar+=j;
    gn++;
  }

  fwrite(&flag,sizeof(flag),1,pfout);
  assert(gn == (iden.ngrupos-1));

  #ifdef STORE_VELOCITIES
  flag = 3*header.npart[1]*sizeof(type_real);
  fwrite(&flag,sizeof(flag),1,pfout);

  npar = gn = 0;
  for(i=1;i<iden.ngrupos;i++)
  {
    j = 0;
    k = Temp.head[i];

    while(k != cp.npart)
    {
      fwrite(&P.vx[k],sizeof(type_real),1,pfout);
      fwrite(&P.vy[k],sizeof(type_real),1,pfout);
      fwrite(&P.vz[k],sizeof(type_real),1,pfout);
      k = Temp.ll[k];
      j++;
    }
    
    assert(j == Temp.npgrup[i]);

    npar+=j;
    gn++;
  }

  fwrite(&flag,sizeof(flag),1,pfout);
  assert(gn == (iden.ngrupos-1));
  #endif    
   
  #ifdef STORE_IDS
  flag = header.npart[1]*sizeof(type_int);
  fwrite(&flag,sizeof(flag),1,pfout);

  npar = gn = 0;
  for(i=1;i<iden.ngrupos;i++)
  {
    j = 0;
    k = Temp.head[i];

    while(k != cp.npart)
    {
      fwrite(&P.id[k],sizeof(type_int),1,pfout);
      k = Temp.ll[k];
      j++;
    }
    
    assert(j == Temp.npgrup[i]);

    npar+=j;
    gn++;
  }

  fwrite(&flag,sizeof(flag),1,pfout);
  assert(gn == (iden.ngrupos-1));
  #endif

  fclose(pfout);

  ///////////////////////////////////////////////////////
  sprintf(filename,"%.2d_%.2f_fof.dat",snap.num,fof[niv]);
  pfout=fopen(filename,"w");
  //////////////////////////////////////////////////////
  for(i=1;i<iden.ngrupos;i++)
  {
    j = 0;
    k = Temp.head[i];
    while(k != cp.npart)
    {
      fprintf(pfout,"%d\n",i);
      k = Temp.ll[k];
      j++;
    }
    assert(j == Temp.npgrup[i]);
  }
  fclose(pfout);
   

  fprintf(stdout,"num de grupos %u num de particulas en grupos %u\n",gn,npar);
  fflush(stdout);

  return;
}

extern void identification(void)
{
  type_int i, j, tid;
  type_int ngrid_old = grid.ngrid;
  type_int *test;  

  gr = (unsigned int **) malloc(nfrac*sizeof(unsigned int *));
  for(j=0;j<nfrac;j++)
    gr[j] = (unsigned int *) malloc(cp.npart*sizeof(unsigned int));

  #ifdef LOCK
    lock = (omp_lock_t *) malloc(cp.npart*sizeof(omp_lock_t));
  #endif
  for(i=0;i<cp.npart;i++)
  {
    #ifdef LOCK
      omp_init_lock(&(lock[i]));
    #endif
    for(j=0; j<nfrac; j++)
    {
      gr[j][i] = i;
    }
  }

  if(nfrac!=1)
  {
    i = (iden.r0[0]>iden.r0[nfrac-1]) ? 0 : nfrac-1; // ORDEN DECRECIENTE O CRECIENTE
  }else{
    i = 0;  
  }

  grid.ngrid = (long)(cp.lbox/iden.r0[i]);
  test = (type_int *) calloc(cp.npart/32 + 1,sizeof(type_int));

  if(grid.ngrid > NGRIDMAX){
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }
  
  grid.nobj = iden.nobj;

  if(ngrid_old != grid.ngrid)
  {
    grid_free();
    grid_init();
    grid_build();
  }

  for(i=0;i<nfrac;i++)
    iden.r0[i] *= iden.r0[i];

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif
    
  fprintf(stdout,"Comienza identificacion.....\n");
  fprintf(stdout,"\nRunning on %d threads\n",NTHREADS);
  fflush(stdout);

  #ifdef LOCK
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,cp,test,lock,stdout)   
  #else
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,cp,test,stdout)   
  #endif
  {
    tid = omp_get_thread_num(); 

    for(i = tid*floor((float)cp.npart/NTHREADS);
    i<(tid==NTHREADS-1 ? cp.npart : (tid+1)*floor((float)cp.npart/NTHREADS));
    i++)
    {
     
      if(i%1000000==0) fprintf(stdout,"%u %u %.4f\n",tid,i,(float)i/(float)cp.npart);

      //#pragma omp taskwait
      //#pragma omp task

      #pragma omp taskgroup
      {

        #ifdef LOCK

          omp_set_lock(&(lock[i]));

          if(TestBit(test,i))
          { 

            omp_unset_lock(&(lock[i]));

          }else{

            SetBit(test,i);
            omp_unset_lock(&(lock[i]));
            busv(i,test);
          } 

        #else

          if(!TestBit(test,i))
          {
            SetBit(test,i);
            busv(i,test);
          }

        #endif
      }
    }

  }  /****** TERMINA SECCION PARALELA ****************************/

  fprintf(stdout,"Sale del paralelo\n"); fflush(stdout);

  #ifdef LOCK
  for(i=0;i<cp.npart;i++) 
    omp_destroy_lock(&(lock[i]));
  free(lock);
  #endif

  for(j=0;j<cp.npart;j++)
  {
    P.x[j] += pmin[0];
    P.y[j] += pmin[1];
    P.z[j] += pmin[2];
  }

  for(j=0;j<nfrac;j++)
  {
    linkedlist(test, gr[j]);
    Write_Groups(j);

    free(Temp.head);
    free(Temp.npgrup);
    free(Temp.ll);
  }
  free(test);

  for(j=0;j<nfrac;j++)
    free(gr[j]);

  fprintf(stdout,"Termino identificacion\n"); fflush(stdout);

}

