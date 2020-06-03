#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "propiedades.h"
#include "variables.h"
#include "allocate.h"
#include "leesnap.h"
#include "grid.h"
#include "iden.h"
#include "io.h"

#define DIV_CEIL(x,y) (x+y-1)/y

static type_int nstep;
static struct iden_st iden;
static struct temporary Temp;
static type_int **gr;
#ifdef LOCK
  static omp_lock_t *lock;
#endif
 
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

static void busv(const type_int ic)
{

  long ixc, iyc, izc, ibox;
  type_int i, niv;
  type_real xx, yy, zz;

#ifdef COLUMN

  ixc  = (long)(P.x[ic]*(type_real)grid.ngrid*(1.f/cp.lbox));
  iyc  = (long)(P.y[ic]*(type_real)grid.ngrid*(1.f/cp.lbox));
  izc  = (long)(P.z[ic]*(type_real)grid.ngrid*(1.f/cp.lbox));

#else

  ixc  = (long)(P[ic].pos[0]*(type_real)grid.ngrid*(1.f/cp.lbox));
  iyc  = (long)(P[ic].pos[1]*(type_real)grid.ngrid*(1.f/cp.lbox));
  izc  = (long)(P[ic].pos[2]*(type_real)grid.ngrid*(1.f/cp.lbox));

#endif

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

        for(i=grid.icell[ibox];i<grid.icell[ibox+1];i++)
        {
          if(ic<i)
          {

#ifdef COLUMN           
            xx = P.x[i] - P.x[ic];
            yy = P.y[i] - P.y[ic];
            zz = P.z[i] - P.z[ic];
#else
            xx = P[i].pos[0] - P[ic].pos[0];
            yy = P[i].pos[1] - P[ic].pos[1];
            zz = P[i].pos[2] - P[ic].pos[2];
#endif
            #ifdef PERIODIC
            xx = ( xx >  cp.lbox*0.5f ) ? xx - cp.lbox : xx ;
            yy = ( yy >  cp.lbox*0.5f ) ? yy - cp.lbox : yy ;
            zz = ( zz >  cp.lbox*0.5f ) ? zz - cp.lbox : zz ;
            xx = ( xx < -cp.lbox*0.5f ) ? xx + cp.lbox : xx ;
            yy = ( yy < -cp.lbox*0.5f ) ? yy + cp.lbox : yy ;
            zz = ( zz < -cp.lbox*0.5f ) ? zz + cp.lbox : zz ;
            #endif

            for(niv=0;niv<nstep;niv++)
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

static void linkedlist(type_int * restrict ar)
{
  type_int i, g;
  
  Temp.ll = (type_int *) calloc(iden.nobj,sizeof(type_int));

  iden.ngrupos = 0;
  for(i=0;i<iden.nobj;i++)
  {
    ar[i] = Raiz(i,ar);
    assert(ar[i]>=i);
    if(Temp.ll[ar[i]] < NPARTMIN)
    {
      Temp.ll[ar[i]]++;
      if(Temp.ll[ar[i]]==NPARTMIN)
      { 
        iden.ngrupos++;
        Temp.ll[ar[i]] = NPARTMIN + iden.ngrupos;
      }
    }
  }

  iden.ngrupos++;  // SUMA UNO;

  Temp.head   = (type_int *) malloc(iden.ngrupos*sizeof(type_int));
  Temp.npgrup = (type_int *) malloc(iden.ngrupos*sizeof(type_int));

  for(i=0;i<iden.ngrupos;i++)
  {
    Temp.head[i]   = iden.nobj;
    Temp.npgrup[i] =  0;
  }

  for(i=0;i<iden.nobj;i++)
  {
    if(Temp.ll[ar[i]]>NPARTMIN)
    {
      ar[i] = Temp.ll[ar[i]] - NPARTMIN;
    }else{
#ifdef COLUMN
      P.sub[i] = ar[i] = 0;
#else
      P[i].sub = ar[i] = 0;
#endif
    }

    g = ar[i];

    #ifdef DEBUG
    assert((g >= 0) && (g < iden.ngrupos));
    #endif
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
    Temp.npgrup[g]++;
  }

  return;
}

static void Write_Groups(const type_int niv)
{
  type_int i,j,k,dim,npar,gn,save_sub;
  type_real dx[3];
  char filename[200];
  FILE *pfout, *pfcentros;
  #ifdef FILE_ASCII
    FILE *pfcentros_ascii;
  #endif
  struct propiedades_st Prop;

  i = iden.ngrupos-1; // LE RESTO UNO POR EL GRUPO 0 PARA ESCRIBIR EN EL ARCHIVO

  ///////////////////////////////////////////////////////
  sprintf(filename,"%.2d_%.2d_level_%.2f_fof.bin",snap.num,niv,fof[niv]);
  pfout=fopen(filename,"w");
  fwrite(&i,sizeof(type_int),1,pfout);
  //////////////////////////////////////////////////////
  sprintf(filename,"%.2d_%.2d_level_%.2f_centros.bin",snap.num,niv,fof[niv]);
  pfcentros=fopen(filename,"w");
  fwrite(&i,sizeof(type_int),1,pfcentros);
  //////////////////////////////////////////////////////
  #ifdef FILE_ASCII
    sprintf(filename,"%.2d_%.2d_level_%.2f_centros.dat",snap.num,niv,fof[niv]);
    pfcentros_ascii=fopen(filename,"w");
  #endif  
  //////////////////////////////////////////////////////

  npar = gn = 0;

  for(i=1;i<iden.ngrupos;i++)
  {

    j = 0;
    k = Temp.head[i];
#ifdef COLUMN           
    save_sub = P.sub[k] == 0 ? i : P.sub[k];
#else
    save_sub = P[k].sub == 0 ? i : P[k].sub;
#endif
    fwrite(&save_sub,sizeof(type_int),1,pfout);
    fwrite(&i,sizeof(type_int),1,pfout);
    fwrite(&Temp.npgrup[i],sizeof(type_int),1,pfout);    

    Prop.npart = Temp.npgrup[i];
    for(dim = 0; dim < 3; dim++)
	  {
		  Prop.pcm[dim] = 0.;
#ifdef STORE_VELOCITIES        
		  Prop.vcm[dim] = 0.;
		  Prop.L[dim]   = 0.;
		  Prop.sig[dim] = 0.;
#endif
	  }
    Prop.pos = (type_real *) malloc(3*Prop.npart*sizeof(type_real));
#ifdef STORE_VELOCITIES        
    Prop.vel = (type_real *) malloc(3*Prop.npart*sizeof(type_real));
#endif

    while(k != iden.nobj)
    {

      // cuidado con el orden {pos[i]-centro} en este caso
#ifdef COLUMN
      Prop.pos[3*j+0] = P.x[k];
      Prop.pos[3*j+1] = P.y[k];
      Prop.pos[3*j+2] = P.z[k];
#ifdef STORE_VELOCITIES        
      Prop.vel[3*j+0] = P.vx[k];
      Prop.vel[3*j+1] = P.vy[k];
      Prop.vel[3*j+2] = P.vz[k];
#endif

      dx[0] = P.x[k] - P.x[Temp.head[i]];
      dx[1] = P.y[k] - P.y[Temp.head[i]];
      dx[2] = P.z[k] - P.z[Temp.head[i]];
#else
      Prop.pos[3*j+0] = P[k].pos[0];
      Prop.pos[3*j+1] = P[k].pos[1];
      Prop.pos[3*j+2] = P[k].pos[2];
#ifdef STORE_VELOCITIES        
      Prop.vel[3*j+0] = P[k].vel[0];
      Prop.vel[3*j+1] = P[k].vel[1];
      Prop.vel[3*j+2] = P[k].vel[2];
#endif

      dx[0] = P[k].pos[0] - P[Temp.head[i]].pos[0];
      dx[1] = P[k].pos[1] - P[Temp.head[i]].pos[1];
      dx[2] = P[k].pos[2] - P[Temp.head[i]].pos[2];
#endif

      for(dim=0; dim<3; dim++)
      {
        #ifdef PERIODIC
        if(dx[dim] >  cp.lbox*0.5)
        {
          dx[dim] -= cp.lbox;
          Prop.pos[3*j+dim] -= cp.lbox;
        }

        if(dx[dim] < -cp.lbox*0.5)
        {
          dx[dim] += cp.lbox;
          Prop.pos[3*j+dim] += cp.lbox;
        }
        #endif

        Prop.pcm[dim] += dx[dim];
#ifdef STORE_VELOCITIES        
        Prop.vcm[dim] += Prop.vel[3*j + dim];
#endif
      }

#ifdef COLUMN
      fwrite(&P.id[k],sizeof(type_int),1,pfout);
      P.sub[k] = i;
#else
      fwrite(&P[k].id,sizeof(type_int),1,pfout);
      P[k].sub = i;
#endif
      k = Temp.ll[k];
      j++;
    }
    
    assert(j == Temp.npgrup[i]);

    for(dim = 0; dim < 3; dim++)
	  {
	  	Prop.pcm[dim] /= (type_real)Prop.npart;
#ifdef STORE_VELOCITIES        
	  	Prop.vcm[dim] /= (type_real)Prop.npart;
#endif

      //	Recenter	
      Prop.pcm[dim] += Prop.pos[dim];

#ifdef CHANGE_POSITION
      Prop.pcm[dim] += pmin[dim];
#endif

#ifdef PERIODIC
      Prop.pcm[dim] = Prop.pcm[dim]<0.0f     ? cp.lbox+Prop.pcm[dim] : Prop.pcm[dim];
      Prop.pcm[dim] = Prop.pcm[dim]>=cp.lbox ? Prop.pcm[dim]-cp.lbox : Prop.pcm[dim];
#endif

    }

    propiedades(&Prop);

    fwrite(&save_sub,sizeof(type_int),1,pfcentros);
    fwrite(&i,sizeof(type_int),1,pfcentros);
    write_properties(pfcentros, Prop);

    #ifdef FILE_ASCII
      fprintf(pfcentros_ascii,"%u %u ",save_sub,i);
      write_properties_ascii(pfcentros_ascii, Prop);
    #endif
     
    free(Prop.pos);
#ifdef STORE_VELOCITIES        
    free(Prop.vel);
#endif

    npar+=j;
    gn++;
  }

  assert(gn == (iden.ngrupos-1));
  fclose(pfout);
  fclose(pfcentros);
  #ifdef FILE_ASCII
    fclose(pfcentros_ascii);
  #endif

  fprintf(stdout,"num de grupos %u num de particulas en grupos %u\n",gn,npar);
  fflush(stdout);

  return;
}

extern void identification(void)
{
  type_int i, j, step, tid;

  for(step=0;step<DIV_CEIL(nfrac,NCUT);step++)
  {
    nstep = step == DIV_CEIL(nfrac,NCUT) - 1 ? nfrac-NCUT*step : NCUT;

    fprintf(stdout,"\n%d Step - %d Levels\n",step,nstep);

    iden.nobj = cp.npart;
    iden.r0 = (double *) malloc(nstep*sizeof(double));
    gr      = (type_int **) malloc(nstep*sizeof(type_int *));

    for(j=0;j<nstep;j++)
    {
      gr[j] = (type_int *) malloc(iden.nobj*sizeof(type_int));
      iden.r0[j]  = fof[NCUT*step+j];
      iden.r0[j] *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0f; //EN KPC

      if(iden.r0[j] <= cp.soft)
      {
        fprintf(stdout,"cambia Linking length = %f \n",iden.r0[j]);
        iden.r0[j] = cp.soft;
      }

      fprintf(stdout,"Linking length %d = %f \n",2*step+j,iden.r0[j]);
    }

    #ifdef LOCK
      lock = (omp_lock_t *) malloc(iden.nobj*sizeof(omp_lock_t));
    #endif
    for(i=0;i<iden.nobj;i++)
    {
#ifdef COLUMN           
      P.sub[i] = step == 0 ? 0 : P.sub[i];
#else
      P[i].sub = step == 0 ? 0 : P[i].sub;
#endif
      #ifdef LOCK
        omp_init_lock(&(lock[i]));
      #endif
      for(j=0; j<nstep; j++)
      {
        gr[j][i] = i;
      }
    }

    grid.ngrid = (long)(cp.lbox/iden.r0[0]);

    if(grid.ngrid > NGRIDMAX)
    {
      fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
      grid.ngrid = NGRIDMAX;
    }
    
    grid.nobj = iden.nobj;
    grid_init();
    grid_build();

    for(j=0;j<nstep;j++)
      iden.r0[j] *= iden.r0[j];

    #ifdef NTHREADS
    omp_set_dynamic(0);
    omp_set_num_threads(NTHREADS);
    #endif
      
    fprintf(stdout,"Comienza identificacion.....\n");
    fprintf(stdout,"Running on %d threads\n",NTHREADS);
    fflush(stdout);

    #ifdef LOCK
      #pragma omp parallel default(none) private(tid,i) \
      shared(P,iden,cp,lock,stdout)   
    #else
      #pragma omp parallel default(none) private(tid,i) \
      shared(P,iden,cp,stdout)   
    #endif
    {
      tid = omp_get_thread_num(); 

      for(i = tid*DIV_CEIL(iden.nobj,NTHREADS);
      i<(tid==NTHREADS-1 ? iden.nobj : (tid+1)*DIV_CEIL(iden.nobj,NTHREADS));
      i++)
      {
       
        if(i%1000000==0) fprintf(stdout,"%u %u %u %.4f\n",tid,i,iden.nobj,(float)i/(float)iden.nobj);

        #pragma omp task
        {
          busv(i);
        }
      }

    }  /****** TERMINA SECCION PARALELA ****************************/

    fprintf(stdout,"Sale del paralelo\n"); fflush(stdout);

    #ifdef LOCK
    for(i=0;i<iden.nobj;i++) 
      omp_destroy_lock(&(lock[i]));
    free(lock);
    #endif

    for(j=0;j<nstep;j++)
    {
      linkedlist(gr[j]);
      Write_Groups(NCUT*step+j);

      free(gr[j]);
      free(Temp.ll);
      free(Temp.head);
      free(Temp.npgrup);
    }

    j = 0;
    for(i=0;i<iden.nobj;i++)
    {
#ifdef COLUMN           
      if(P.sub[i] != 0)
      {
        P.x[j] = P.x[i];
        P.y[j] = P.y[i];
        P.z[j] = P.z[i];
        #ifdef STORE_VELOCITIES
          P.vx[j] = P.vx[i];
          P.vy[j] = P.vy[i];
          P.vz[j] = P.vz[i];
        #endif
        #ifdef STORE_IDS
          P.id[j] = P.id[i];
        #endif
        P.sub[j] = P.sub[i];
        j++;
      }
#else
      if(P[i].sub != 0)
      {
        P[j] = P[i];
        j++;
      }
#endif
  }
  
    cp.npart = j;
    if(!reallocate_particles(&P, cp.npart))  exit(1);

    free(iden.r0);
    grid_free();
  }

  return;
}

