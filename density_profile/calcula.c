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

#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define DIV_CEIL(i,j) (i+j-1)/j

//////////////////////////////
static type_real *rcil2; 
static type_real r_aux; // RADIO DE BUSQUEDA EN Kpc
//////////////////////////////

static void set_name(char * name, const int id)
{
  struct stat st = {0};

  #ifdef ROTACION
    sprintf(name,"./%.2d_rotando_excluyendo_%.2f_rvir/",TYPE_FLAG,RVIR_FACTOR);
  #else
    #ifdef ESFERA
      sprintf(name,"./%.2d_excluyendo_%.2f_rvir/".TYPE_FLAG,RVIR_FACTOR);
    #else
      sprintf(name,"./%.2d_sin_esfera_excluyendo_%.2f_rvir/",TYPE_FLAG,RVIR_FACTOR);
    #endif
  #endif

  if(stat(name, &st) == -1) 
  {
    mkdir(name, 0700);
  }

  sprintf(name,"%s%.2d_type_%.2d",name,id,TYPE_FLAG);

  #ifdef BIN_LOG
    sprintf(name,"%s_%.2d_%.2f_Mpc_LOG",name,ncil,pow(10.0,RAUX)/1000.);
  #else
    sprintf(name,"%s_%.2d_%.2f_Mpc_LIN",name,ncil,RAUX/1000.);
  #endif

  #ifdef SIN_REPETICION 
    sprintf(name,"%s_sin_repeticion",name);
  #else
    sprintf(name,"%s_con_repeticion",name);
  #endif

  #ifdef BINARIO
    sprintf(name,"%s.bin",name);
  #else
    sprintf(name,"%s.dat",name);
  #endif
  

  return;
}

static void set_name_prop(char * name, const int id)
{
  struct stat st = {0};

  #ifdef ROTACION
    sprintf(name,"./%.2d_rotando_excluyendo_%.2f_rvir/",TYPE_FLAG,RVIR_FACTOR);
  #else
    #ifdef ESFERA
      sprintf(name,"./%.2d_excluyendo_%.2f_rvir/",TYPE_FLAG,RVIR_FACTOR);
    #else
      sprintf(name,"./%.2d_sin_esfera_excluyendo_%.2f_rvir/",TYPE_FLAG,RVIR_FACTOR);
    #endif
  #endif

  if(stat(name, &st) == -1) 
  {
    mkdir(name, 0700);
  }

  sprintf(name,"%s%.2d_type_%.2d_rlim_%.2f_Mpc",name,id,TYPE_FLAG,RSPH/1000.);

  #ifdef SIN_REPETICION 
    sprintf(name,"%s_propiedades_sin_repeticion",name);
  #else
    sprintf(name,"%s_propiedades_con_repeticion",name);
  #endif

  #ifdef BINARIO
    sprintf(name,"%s.bin",name);
  #else
    sprintf(name,"%s.dat",name);
  #endif

  return;
}

#ifdef BIN_LOG

  static void logspace(type_real *rcil2, const type_real max, const type_real min, const type_int bins) 
  {
  
    type_int i;
    type_real end   = max;
    type_real start = min;
    type_real delta = (end - start) / (bins-1);
  
    for(i=0;i<bins; i++)
    {
      rcil2[i] = pow(10.0,(start + delta * i));
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
  
    for(i=0;i<bins; i++)
    {
      rcil2[i] = start + delta * i;
      rcil2[i] *= rcil2[i];
    }  
  
    return;
  }

#endif

static int check_exists(char * filename)
{
  FILE *pf;
  // try to open file to read
  
  if(!(pf=fopen(filename,"r")))
  {
    sprintf(message,"not exist file %s\n",filename);RED(message);
    fflush(stdout);
    return 0;
  }else{
    sprintf(message,"exist file %s\n",filename);GREEN(message);
    fflush(stdout);
    fclose(pf);
    return 1;
  } 
}

extern void limpia_calculados()
{
  char filename[200];
  type_int i,k,c,l;

  c = l = 0;
  for(i=0;i<cp.nseg;i++)
  {
    set_name(filename,Seg[i].id);

    if(check_exists(filename) == 1) continue;

    for(k=Seg[i].start;k<Seg[i].start+Seg[i].size;k++)
    {
      Gr[c].Pos[0] = Gr[k].Pos[0];
      Gr[c].Pos[1] = Gr[k].Pos[1];
      Gr[c].Pos[2] = Gr[k].Pos[2];
      
      c++;
    }

    Seg[i].start = c-Seg[i].size;
    Seg[l] = Seg[i];

    l++;
  }

  cp.nseg = l;
  cp.ngrup = c;

  Seg = (struct segmentstd *) realloc(Seg,cp.nseg*sizeof(struct segmentstd));
  Gr  = (struct grup_data  *) realloc(Gr,cp.ngrup*sizeof(struct grup_data));

  fprintf(stdout,"Filamentos CUT  %d\n",cp.nseg);
  fprintf(stdout,"Nodos      CUT  %d\n",cp.ngrup);
  fprintf(stdout,"Segmentos  CUT  %d\n",cp.ngrup-cp.nseg);
  fflush(stdout);

  return;
}

#ifdef SIN_REPETICION

  static void setRandomSeed(unsigned short *seed16)
  {
    struct timeval tv;			// Time of day, for getting usecs
    struct timezone tz;			// Dummy for gettimeofday
    gettimeofday (&tv, &tz);
  
    double wideSeed = (unsigned long)tv.tv_sec * 1000000.0 + tv.tv_usec;
    const double divisor = 1ul << 16;
    unsigned long topBits = (unsigned long)(wideSeed / divisor / divisor);
    seed16[2] = (unsigned short)topBits;
    wideSeed -= topBits * divisor * divisor;
    seed16[1] = (unsigned short)(wideSeed / divisor);
    wideSeed -= seed16[1] * divisor;
    seed16[0] = (unsigned short)wideSeed;
    
    return;
  }
  
  static void decode(const long ijk, const long ngrid, long *ix, long *iy, long *iz)
  {
    long ijkt;
  
  	*ix   = ijk/(ngrid*ngrid);	
  	ijkt =  ijk - ngrid*ngrid*(*ix);
  	*iy   = ijkt/ngrid;
  	*iz   = ijkt-(*iy)*ngrid;
  
    return; 
  }
  
  static int compare_ids(const void *a, const void *b)
  { 
    int *c, *d;
    
    c = (int *) a;
    d = (int *) b;
  
    return (*c - *d);
  }

#ifndef SERIAL
  
  static void calc_profile_sin_repeticion(void)
  {
    int i, j, k, id_min;
    long ixx, iyy, izz;
    long ixc, iyc, izc, naux;
    type_real dis, rsep, rsep_min, dis_min, rho;
    type_real vaux[3], delta[3], Pos_cent[3];
    FILE *pf;
    long  *work_data;
    type_int Tid, cont;
    type_int **dens_rand;
    type_int **dens_data;
    type_real *tmp_matrix_mean_profile;
    type_real *tmp_matrix_mean_cil;
    type_int  mass;
    type_real mu;
    unsigned short **seed16;
    char filename[200];
    const type_int NTHREADS = omp_get_max_threads();
  
    seed16 = (unsigned short **)  malloc(NTHREADS*sizeof(unsigned short *));
    for(k=0;k<NTHREADS;k++)
      seed16[k] = (unsigned short *) malloc(3*sizeof(unsigned short));
  
    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) \
    private(i,j,k,id_min,rsep,dis,dis_min,rsep_min,\
    vaux,delta,Pos_cent,Tid,rho,ixc,iyc,izc,ixx,iyy,izz,\
    dens_data,dens_rand,work_data,tmp_matrix_mean_cil,tmp_matrix_mean_profile,mu,mass,\
    cont,naux,pf,filename) shared(cp,ncil,RAUX,RINIT,r_aux,\
    rcil2,Seg,Gr,P,grid,seed16,stdout)
    for(i=0;i<cp.nseg;i++)
    {
  
      Tid = omp_get_thread_num();
  
      set_name(filename,Seg[i].id);
  
      if(check_exists(filename) == 1) continue;
  
      // INIT ARRAYS //////////////////////////////////

      pf=fopen("/dev/urandom","r");
      for(k=0;k<3;k++)
        fread(&seed16[Tid][k],sizeof(unsigned short),1,pf);
      fclose(pf);     

      work_data = (long *)  malloc(27*(Seg[i].size)*sizeof(long));
      dens_rand = (type_int **) calloc((Seg[i].size-1),sizeof(type_int *));
      dens_data = (type_int **) calloc((Seg[i].size-1),sizeof(type_int *));
      tmp_matrix_mean_profile = (type_real  *) calloc(ncil,sizeof(type_real));
      tmp_matrix_mean_cil     = (type_real  *) calloc(2*ncil,sizeof(type_real));
  
      for(k=0;k<Seg[i].size-1;k++)
      {
        dens_rand[k] = (type_int *) calloc(ncil,sizeof(type_int));
        dens_data[k] = (type_int *) calloc(ncil,sizeof(type_int));
      }

      mu = 0.0;
      mass = 0.0;

      ////////////////////////////////////////////////
  
      fprintf(stdout,"size %d seg %d/%d fil %d/%d\n",Seg[i].size-1,Seg[i].start-i,cp.ngrup-cp.nseg,i,cp.nseg);
      fflush(stdout);
  
      cont = 0;

      //////////////////////////////////////////////
      ixc  = (long)((Gr[Seg[i].start].Pos[0])*(type_real)grid.ngrid*(1.f/cp.lbox));
      iyc  = (long)((Gr[Seg[i].start].Pos[1])*(type_real)grid.ngrid*(1.f/cp.lbox));
      izc  = (long)((Gr[Seg[i].start].Pos[2])*(type_real)grid.ngrid*(1.f/cp.lbox));
  
      #ifndef PERIODIC
        for(ixx = ((ixc-1<0) ? 0 : ixc-1); ixx <= ((ixc+1 >= grid.ngrid) ? grid.ngrid-1 : ixc+1); ixx++)
      #else
        for(ixx = ixc-1; ixx <= ixc+1; ixx++)
      #endif
      {
        #ifndef PERIODIC
          for(iyy = ((iyc-1<0) ? 0 : iyc-1); iyy <= ((iyc+1 >= grid.ngrid) ? grid.ngrid-1 : iyc+1); iyy++)
        #else
          for(iyy = iyc-1 ; iyy <= iyc+1 ; iyy++)
        #endif
        {
          #ifndef PERIODIC
            for(izz = ((izc-1<0) ? 0 : izc-1); izz <= ((izc+1 >= grid.ngrid) ? grid.ngrid-1 : izc+1); izz++)
          #else
            for(izz = izc-1 ; izz <= izc+1 ; izz++)
          #endif
          {
  
          	#ifdef PERIODIC
            	work_data[cont++] = ( ( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) )   * grid.ngrid +\
  	  			                        ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ) ) * grid.ngrid +\
                                    ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) );
            #else
        			work_data[cont++] = (ixx * grid.ngrid + iyy) * grid.ngrid + izz;
            #endif
         } // izz
        } // iyy
      } // ixx         

      dis = 0.0;
      for(k=Seg[i].start+1;k<Seg[i].start+Seg[i].size;k++)
      {

        delta[0] = Gr[k].Pos[0] - Gr[k-1].Pos[0];
        delta[1] = Gr[k].Pos[1] - Gr[k-1].Pos[1];
        delta[2] = Gr[k].Pos[2] - Gr[k-1].Pos[2];
  
        #ifdef PERIODIC
        if(delta[0]> cp.lbox) delta[0] -= cp.lbox;
        if(delta[0]< 0.0)     delta[0] += cp.lbox;
        if(delta[1]> cp.lbox) delta[1] -= cp.lbox;
        if(delta[1]< 0.0)     delta[1] += cp.lbox;
        if(delta[2]> cp.lbox) delta[2] -= cp.lbox;
        if(delta[2]< 0.0)     delta[2] += cp.lbox;
        #endif
        
        dis += sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);

        ///////////////////////////////////////////////

        Pos_cent[0] = Gr[k].Pos[0];
        Pos_cent[1] = Gr[k].Pos[1];
        Pos_cent[2] = Gr[k].Pos[2];
  
        #ifdef PERIODIC
        if(Pos_cent[0]> cp.lbox) Pos_cent[0] -= cp.lbox;
        if(Pos_cent[0]< 0.0)     Pos_cent[0] += cp.lbox;
        if(Pos_cent[1]> cp.lbox) Pos_cent[1] -= cp.lbox;
        if(Pos_cent[1]< 0.0)     Pos_cent[1] += cp.lbox;
        if(Pos_cent[2]> cp.lbox) Pos_cent[2] -= cp.lbox;
        if(Pos_cent[2]< 0.0)     Pos_cent[2] += cp.lbox;
        #endif
  
      	ixc  = (long)((Pos_cent[0])*(type_real)grid.ngrid*(1.f/cp.lbox));
      	iyc  = (long)((Pos_cent[1])*(type_real)grid.ngrid*(1.f/cp.lbox));
      	izc  = (long)((Pos_cent[2])*(type_real)grid.ngrid*(1.f/cp.lbox));
  
        #ifndef PERIODIC
          for(ixx = ((ixc-1<0) ? 0 : ixc-1); ixx <= ((ixc+1 >= grid.ngrid) ? grid.ngrid-1 : ixc+1); ixx++)
        #else
          for(ixx = ixc-1; ixx <= ixc+1; ixx++)
        #endif
        {
          #ifndef PERIODIC
            for(iyy = ((iyc-1<0) ? 0 : iyc-1); iyy <= ((iyc+1 >= grid.ngrid) ? grid.ngrid-1 : iyc+1); iyy++)
          #else
            for(iyy = iyc-1 ; iyy <= iyc+1 ; iyy++)
          #endif
          {
            #ifndef PERIODIC
              for(izz = ((izc-1<0) ? 0 : izc-1); izz <= ((izc+1 >= grid.ngrid) ? grid.ngrid-1 : izc+1); izz++)
            #else
              for(izz = izc-1 ; izz <= izc+1 ; izz++)
            #endif
            {
  
            	#ifdef PERIODIC
              	work_data[cont++] = ( ( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) )   * grid.ngrid +\
  	    			                        ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ) ) * grid.ngrid +\
                                      ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) );
              #else
          			work_data[cont++] = (ixx * grid.ngrid + iyy) * grid.ngrid + izz;
              #endif
           } // izz
          } // iyy
        } // ixx         
      } // cierra el for de k
  
      qsort(work_data,cont,sizeof(long),&compare_ids);
  
      ixx = 0;   
      for(iyy=0; iyy<cont-1; iyy++) 
        if(work_data[iyy] != work_data[iyy+1]) 
          work_data[ixx++] = work_data[iyy]; 
  
      work_data[ixx++] = work_data[cont-1];   
      cont = ixx;

      /// DISTANCIA ENTRE LAS PUNTAS      
      delta[0] = Gr[Seg[i].start+Seg[i].size-1].Pos[0]-Gr[Seg[i].start].Pos[0];
      delta[1] = Gr[Seg[i].start+Seg[i].size-1].Pos[1]-Gr[Seg[i].start].Pos[1];
      delta[2] = Gr[Seg[i].start+Seg[i].size-1].Pos[2]-Gr[Seg[i].start].Pos[2];
  
      #ifdef PERIODIC
      if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
      if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
      if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
      if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
      if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
      if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
      #endif
    
      // distancia entre los nodos
      rsep = sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]); 

      //ELONGACION
      rsep /= dis; // DIVIDO POR LA LONGITUD TOTAL

      dis /= (type_real)Seg[i].size-1;   // LONGITUD PROMEDIO DE SEGMENTOS

      dis *= rsep;                       // LA IDEA ES QUE MIENTRAS MAS ELONGADO
                                         // VOY A TIRAR MAS RANDOM

      dis *= M_PI*(rcil2[1] - rcil2[0]); // VOLUMEN DEL PRIMER CILINDRO
      dis = pow(r_aux,3.0)/dis;          // VOLUMEN DE LA CELDA    
      izz = NLIMITE*(long)dis;           // CANTIDAD DE RANDOM POR CELDA
      naux = izz>(NRANDOM*(long)((cp.npart + cont-1)/cont)) ? \
      (NRANDOM*(long)((cp.npart + cont-1)/cont)) : izz;   // CANTIDAD DE RANDOM si se pasa de un valor critico lo cambio

      for(iyy=0; iyy<cont; iyy++) 
      {
        ixx = grid.llirst[work_data[iyy]];
  
        while(ixx != grid.nobj)
        {
           id_min  = -1;
           dis_min = rsep_min = 1e26;
  
           Pos_cent[0] = P[ixx].Pos[0];
           Pos_cent[1] = P[ixx].Pos[1];
           Pos_cent[2] = P[ixx].Pos[2];
  
           // EL ULTIMO NODO
           delta[0] = Pos_cent[0] - Gr[Seg[i].start+Seg[i].size-1].Pos[0];
           delta[1] = Pos_cent[1] - Gr[Seg[i].start+Seg[i].size-1].Pos[1];
           delta[2] = Pos_cent[2] - Gr[Seg[i].start+Seg[i].size-1].Pos[2];
  
           #ifdef PERIODIC
           if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
           if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
           if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
           if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
           if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
           if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
           #endif
  
           dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
           if(dis < Seg[i].Rvir_2[1]) // uso el cuadrado  
           {
             ixx = grid.ll[ixx];  
             continue;
           }
           
           // EL PRIMER NODO
           delta[0] = Pos_cent[0] - Gr[Seg[i].start].Pos[0];
           delta[1] = Pos_cent[1] - Gr[Seg[i].start].Pos[1];
           delta[2] = Pos_cent[2] - Gr[Seg[i].start].Pos[2];
  
           #ifdef PERIODIC
           if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
           if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
           if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
           if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
           if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
           if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
           #endif
  
           dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
           if(dis < Seg[i].Rvir_2[0]) // uso el cuadrado  
           {
             ixx = grid.ll[ixx];  
             continue;
           }

           #ifdef ESFERA
           type_int in_sph_ext = 0;
           #endif
 
           vaux[0] = Gr[Seg[i].start+1].Pos[0]-Gr[Seg[i].start].Pos[0];
           vaux[1] = Gr[Seg[i].start+1].Pos[1]-Gr[Seg[i].start].Pos[1];
           vaux[2] = Gr[Seg[i].start+1].Pos[2]-Gr[Seg[i].start].Pos[2];
  
           #ifdef PERIODIC
           if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
           if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
           if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
           if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
           if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
           if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
           #endif
           rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
           vaux[0] /= rsep;
           vaux[1] /= rsep;
           vaux[2] /= rsep;
  
           // Reutilizo DELTA de lo calculado
           dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];

           // SI ENTRA EN MI CILINDRO
           if(dis>0.0f && dis<rsep) 
           {
             dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;
  
             if(dis>0.0)
             {
               dis = sqrt(dis);
  
               #ifdef BIN_LOG
               dis = log10(dis);
               #endif
  
               if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
               { 
             	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                 if(j < ncil && j >= 0)
                 {
                   id_min = j; // queda igual sumo en el k = 0
                   dis_min = dis;
                   rsep_min = rsep;
                 }
               }
             }
           }

           #ifdef ESFERA
           if(dis>rsep) in_sph_ext = 1;
           #endif

           for(k=Seg[i].start+2;k<Seg[i].start+Seg[i].size;k++)
           {  
             vaux[0] = Gr[k].Pos[0]-Gr[k-1].Pos[0];
             vaux[1] = Gr[k].Pos[1]-Gr[k-1].Pos[1];
             vaux[2] = Gr[k].Pos[2]-Gr[k-1].Pos[2];
  
             #ifdef PERIODIC
             if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
             if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
             if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
             if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
             if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
             if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
             #endif
             rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
             vaux[0] /= rsep;
             vaux[1] /= rsep;
             vaux[2] /= rsep;
  
             delta[0] = Pos_cent[0] - Gr[k-1].Pos[0];
             delta[1] = Pos_cent[1] - Gr[k-1].Pos[1];
             delta[2] = Pos_cent[2] - Gr[k-1].Pos[2];
  
             #ifdef PERIODIC
             if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
             if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
             if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
             if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
             if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
             if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
             #endif
           
             dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];
  
             // SI ENTRA EN MI CILINDRO
             if(dis>0.0f && dis<rsep) 
             {
               dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;
               
               if(dis>0.0)
               {
                 dis = sqrt(dis);
  
                 #ifdef BIN_LOG
                 dis = log10(dis);
                 #endif
  
                 if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                 { 
               	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                   if(j < ncil && j >= 0)
                   {
                     id_min = (k-1-Seg[i].start)*ncil+j;
                     dis_min = dis;                     
                     rsep_min = rsep;
                   }
                 }
               }
  
            #ifdef ESFERA

             in_sph_ext = 0;

             }else if(dis<0.0f){

               if(in_sph_ext==1) 
               {
                 dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
  
                 if(dis>0.0)
                 {
                   dis = sqrt(dis);
  
                   #ifdef BIN_LOG
                   dis = log10(dis);
                   #endif
  
                   if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                   { 
               	    j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                     if(j < ncil && j >= 0)
                     {
                       id_min = (k-1-Seg[i].start)*ncil+j;
                       dis_min = dis;
                       rsep_min = rsep;
                     }
                   }
                 }                           
               }

               in_sph_ext = 0;

             }else{

               in_sph_ext = 1;

             #endif
             }

           } // cierra el for de k
  
           if(id_min != -1)
           {
             k = id_min/ncil;
             j = id_min - ncil*k;    
             dens_data[k][j] += 1;

             if(dis_min <= RSPH)
             {
               mu   += 1.0/rsep_min;
               mass += 1;
             }
           }
           
           ixx = grid.ll[ixx];  
  
         } //fin lazo particulas del grid
        
      } //fin lazo sobre celdas
  
//////////////naux = NRANDOM*((naux + cont-1)/cont);
  
      for(iyy=0;iyy<cont;iyy++)
      {
  
        seed48(seed16[Tid]); // setea la semilla   
  
        decode(work_data[iyy],grid.ngrid,&ixc,&iyc,&izc); 
  
        for(ixx=0;ixx<naux;ixx++)
        {

          id_min  = -1;
          dis_min = 1e26;

          Pos_cent[0] = (drand48()+(type_real)ixc)*r_aux;
          Pos_cent[1] = (drand48()+(type_real)iyc)*r_aux;
          Pos_cent[2] = (drand48()+(type_real)izc)*r_aux;  
 
          #ifdef PERIODIC
          if(Pos_cent[0]> cp.lbox) Pos_cent[0] -= cp.lbox;
          if(Pos_cent[0]< 0.0)     Pos_cent[0] += cp.lbox;
          if(Pos_cent[1]> cp.lbox) Pos_cent[1] -= cp.lbox;
          if(Pos_cent[1]< 0.0)     Pos_cent[1] += cp.lbox;
          if(Pos_cent[2]> cp.lbox) Pos_cent[2] -= cp.lbox;
          if(Pos_cent[2]< 0.0)     Pos_cent[2] += cp.lbox;
          #endif
  
          // EL ULTIMO NODO
          delta[0] = Pos_cent[0] - Gr[Seg[i].start+Seg[i].size-1].Pos[0];
          delta[1] = Pos_cent[1] - Gr[Seg[i].start+Seg[i].size-1].Pos[1];
          delta[2] = Pos_cent[2] - Gr[Seg[i].start+Seg[i].size-1].Pos[2];
  
          #ifdef PERIODIC
          if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
          if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
          if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
          if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
          if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
          if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
          #endif
  
          dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
          if(dis < Seg[i].Rvir_2[1]) // uso el cuadrado  
            continue;
          
          // EL PRIMER NODO
          delta[0] = Pos_cent[0] - Gr[Seg[i].start].Pos[0];
          delta[1] = Pos_cent[1] - Gr[Seg[i].start].Pos[1];
          delta[2] = Pos_cent[2] - Gr[Seg[i].start].Pos[2];
  
          #ifdef PERIODIC
          if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
          if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
          if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
          if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
          if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
          if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
          #endif
  
          dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
          if(dis < Seg[i].Rvir_2[0]) // uso el cuadrado  
            continue;
  
          vaux[0] = Gr[Seg[i].start+1].Pos[0]-Gr[Seg[i].start].Pos[0];
          vaux[1] = Gr[Seg[i].start+1].Pos[1]-Gr[Seg[i].start].Pos[1];
          vaux[2] = Gr[Seg[i].start+1].Pos[2]-Gr[Seg[i].start].Pos[2];
  
          #ifdef PERIODIC
          if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
          if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
          if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
          if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
          if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
          if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
          #endif
          rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
          vaux[0] /= rsep;
          vaux[1] /= rsep;
          vaux[2] /= rsep;

          #ifdef ESFERA
          type_int in_sph_ext = 0;
          #endif

          // Reutilizo DELTA de lo calculado
          dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];
  
          // SI ENTRA EN MI CILINDRO
          if(dis>0.0f && dis<rsep) 
          {
            dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;
  
            if(dis>0.0)
            {
              dis = sqrt(dis);
  
              #ifdef BIN_LOG
              dis = log10(dis);
              #endif
  
              if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
              { 
            	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
  
                if(j < ncil && j >= 0)
                {
                  id_min = j; // queda igual sumo en el k = 0
                  dis_min = dis;
                }
              }
            }
          }

          #ifdef ESFERA
          if(dis>rsep) in_sph_ext = 1;
          #endif

          for(k=Seg[i].start+2;k<Seg[i].start+Seg[i].size;k++)
          {  
            vaux[0] = Gr[k].Pos[0]-Gr[k-1].Pos[0];
            vaux[1] = Gr[k].Pos[1]-Gr[k-1].Pos[1];
            vaux[2] = Gr[k].Pos[2]-Gr[k-1].Pos[2];
  
            #ifdef PERIODIC
            if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
            if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
            if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
            if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
            if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
            if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
            #endif
            rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
            vaux[0] /= rsep;
            vaux[1] /= rsep;
            vaux[2] /= rsep;
  
            delta[0] = Pos_cent[0] - Gr[k-1].Pos[0];
            delta[1] = Pos_cent[1] - Gr[k-1].Pos[1];
            delta[2] = Pos_cent[2] - Gr[k-1].Pos[2];
  
            #ifdef PERIODIC
            if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
            if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
            if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
            if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
            if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
            if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
            #endif
          
            dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];
  
            // SI ENTRA EN MI CILINDRO
            if(dis>0.0f && dis<rsep) 
            {
              dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;

              if(dis>0.0)
              {
                dis = sqrt(dis);
  
                #ifdef BIN_LOG
                dis = log10(dis);
                #endif
  
                if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                { 
              	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
  
                  if(j < ncil && j >= 0)
                  {
                    id_min = (k-1-Seg[i].start)*ncil+j;
                    dis_min = dis;
                  }
                }
              }
            
            #ifdef ESFERA
             in_sph_ext = 0;

            }else if(dis<0.0f){

              if(in_sph_ext==1) 
              {
                dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
  
                if(dis>0.0)
                {
                  dis = sqrt(dis);
  
                  #ifdef BIN_LOG
                  dis = log10(dis);
                  #endif
  
                  if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                  { 
              	    j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                    if(j < ncil && j >= 0)
                    {
                      id_min = (k-1-Seg[i].start)*ncil+j;
                      dis_min = dis;
                    }
                  }
                }                           
              }                           

              in_sph_ext = 0;

             }else{

               in_sph_ext = 1;

              #endif
            }
  
          } // cierra el for de k
  
          if(id_min != -1)
          {
            k = id_min/ncil;
            j = id_min - ncil*k;    
            dens_rand[k][j] += 1;
          }
  
        } // cierra el loop sobre las celdas
  
        setRandomSeed(seed16[Tid]);
  
      } // cierra la cantidad de random      
      
      naux *= (long)cont;                // ntotal
      rho   = (type_real)cont*pow(r_aux,3.0); // vol tot
      rho   = (type_real)naux/rho;            // dens tot
      //rho  /= (type_real)(Seg[i].size-1);     // normalizacion por segmento
  
      cont = 0;
      for(k=0;k<Seg[i].size-1;k++)
      {
        id_min = 1;
        for(j=0;j<ncil;j++)
        {
          id_min *= (dens_rand[k][j]>0 ? 1 : 0);
          tmp_matrix_mean_cil[j]      += (type_real)dens_data[k][j];
          tmp_matrix_mean_cil[ncil+j] += (type_real)dens_rand[k][j];
        }
  
        if(id_min == 1)
        {
          for(j=0;j<ncil;j++)
            tmp_matrix_mean_profile[j] += rho*((type_real)dens_data[k][j]/(type_real)dens_rand[k][j]);
          cont++;
        }
  
        free(dens_rand[k]);
        free(dens_data[k]);
      }
  
      pf = fopen(filename,"w");
  
      #ifdef BINARIO
        fwrite(&ncil,sizeof(type_int),1,pf);        
      #endif
  
      for(j=0;j<ncil;j++)
      {
        dis = 0.5*(sqrt(rcil2[j])+sqrt(rcil2[j+1]))/1000.0; // En Mpc
        tmp_matrix_mean_profile[j] = tmp_matrix_mean_profile[j]*pow(1000.0,3);           // En Mpc
      	tmp_matrix_mean_profile[j] /= (type_real)(cont);     // normalizacion por segmento
  
        tmp_matrix_mean_cil[j]  = rho*(tmp_matrix_mean_cil[j]/tmp_matrix_mean_cil[ncil+j]);
        tmp_matrix_mean_cil[j] *= pow(1000.0,3);           // En Mpc
        // tmp_matrix_mean_cil[j] /= (type_real)(Seg[i].size-1);
  
        #ifdef BINARIO
          fwrite(&dis,sizeof(type_real),1,pf);        
          fwrite(&tmp_matrix_mean_profile[j],sizeof(type_real),1,pf);        
        #else
          fprintf(pf,"%d %d %f %f %f\n",cont,Seg[i].size-1,dis,tmp_matrix_mean_profile[j],tmp_matrix_mean_cil[j]);
        #endif
  
      }
      fclose(pf);
 
      set_name_prop(filename,Seg[i].id);
      pf = fopen(filename,"w");

      mu  *= (cp.Mpart*1000.0);        // En Mass[10^10 Msol / h]/Mpc
      rsep = cp.Mpart*(type_real)mass; // En Mass[10^10 Msol / h]
      dis  = mu/(pow(RSPH/1000.0,2));  // En Mass[10^10 Msol / h]/Mpc^3

      #ifdef BINARIO
        fwrite(&rsep,sizeof(type_real),1,pf);        
        fwrite(&mu,sizeof(type_real),1,pf);        
        fwrite(&dis,sizeof(type_real),1,pf);        
      #else
        fprintf(pf,"%f %f %f\n",rsep,mu,dis);
      #endif
      fclose(pf);

      free(work_data);
      free(dens_rand);
      free(dens_data);
      free(tmp_matrix_mean_cil);
      free(tmp_matrix_mean_profile);
  
    }// cierro el for sobre los filamentos
  
  ////////////////////////////////////////////////////////////////////////////////////////////  
  
    return;
  
  }

#else

  static void calc_profile_sin_repeticion(void)
  {
    int i, j, k, id_min;
    long ixx, iyy, izz;
    long ixc, iyc, izc, naux;
    type_real dis, rsep, rsep_min, dis_min, rho;
    type_real vaux[3], delta[3], Pos_cent[3];
    FILE *pf;
    long  *work_data;
    type_int Tid, cont;
    type_int **dens_rand;
    type_int **dens_data;
    type_real *tmp_matrix_mean_profile;
    type_real *tmp_matrix_mean_cil;
    type_real *Pos_rand;
    type_int  mass;
    type_real mu;
    unsigned short *seed16;
    char filename[200];
    const type_int NTHREADS = omp_get_max_threads();
  
    seed16 = (unsigned short *) malloc(3*sizeof(unsigned short));  
    Pos_rand = (type_real *) malloc((3*NRANDOM*(long unsigned)cp.npart)*sizeof(type_real));

    pf=fopen("/dev/urandom","r");
    for(k=0;k<3;k++)
      fread(&seed16[k],sizeof(unsigned short),1,pf);
    fclose(pf);     

    seed48(seed16); // setea la semilla   

    naux = (NRANDOM*(long)cp.npart);
    for(izz=0;izz<naux;izz++)
      for(k=0;k<3;k++)
        Pos_rand[3*izz+(long)k] = drand48();

    exit(0);

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) \
    private(i,j,k,id_min,rsep,dis,dis_min,rsep_min,\
    vaux,delta,Pos_cent,Tid,rho,ixc,iyc,izc,ixx,iyy,izz,\
    dens_data,dens_rand,work_data,tmp_matrix_mean_cil,tmp_matrix_mean_profile,mu,mass,\
    cont,naux,pf,filename) shared(cp,ncil,RAUX,RINIT,r_aux,\
    rcil2,Seg,Gr,P,grid,Pos_rand,stdout)
    for(i=0;i<cp.nseg;i++)
    {
  
      set_name(filename,Seg[i].id);

      if(check_exists(filename) == 1) continue;
  
      // INIT ARRAYS //////////////////////////////////

      work_data = (long *)  malloc(27*(Seg[i].size)*sizeof(long));
      dens_rand = (type_int **) calloc((Seg[i].size-1),sizeof(type_int *));
      dens_data = (type_int **) calloc((Seg[i].size-1),sizeof(type_int *));
      tmp_matrix_mean_profile = (type_real  *) calloc(ncil,sizeof(type_real));
      tmp_matrix_mean_cil     = (type_real  *) calloc(2*ncil,sizeof(type_real));
  
      for(k=0;k<Seg[i].size-1;k++)
      {
        dens_rand[k] = (type_int *) calloc((NTHREADS*ncil),sizeof(type_int));
        dens_data[k] = (type_int *) calloc((NTHREADS*ncil),sizeof(type_int));
      }

      mu = 0.0;
      mass = 0.0;

      ////////////////////////////////////////////////
  
      fprintf(stdout,"size %d seg %d/%d fil %d/%d\n",Seg[i].size-1,Seg[i].start-i,cp.ngrup-cp.nseg,i,cp.nseg);
      fflush(stdout);
  
      cont = 0;

      //////////////////////////////////////////////
      ixc  = (long)((Gr[Seg[i].start].Pos[0])*(type_real)grid.ngrid*(1.f/cp.lbox));
      iyc  = (long)((Gr[Seg[i].start].Pos[1])*(type_real)grid.ngrid*(1.f/cp.lbox));
      izc  = (long)((Gr[Seg[i].start].Pos[2])*(type_real)grid.ngrid*(1.f/cp.lbox));
  
      #ifndef PERIODIC
        for(ixx = ((ixc-1<0) ? 0 : ixc-1); ixx <= ((ixc+1 >= grid.ngrid) ? grid.ngrid-1 : ixc+1); ixx++)
      #else
        for(ixx = ixc-1; ixx <= ixc+1; ixx++)
      #endif
      {
        #ifndef PERIODIC
          for(iyy = ((iyc-1<0) ? 0 : iyc-1); iyy <= ((iyc+1 >= grid.ngrid) ? grid.ngrid-1 : iyc+1); iyy++)
        #else
          for(iyy = iyc-1 ; iyy <= iyc+1 ; iyy++)
        #endif
        {
          #ifndef PERIODIC
            for(izz = ((izc-1<0) ? 0 : izc-1); izz <= ((izc+1 >= grid.ngrid) ? grid.ngrid-1 : izc+1); izz++)
          #else
            for(izz = izc-1 ; izz <= izc+1 ; izz++)
          #endif
          {
  
          	#ifdef PERIODIC
            	work_data[cont++] = ( ( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) )   * grid.ngrid +\
  	  			                        ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ) ) * grid.ngrid +\
                                    ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) );
            #else
        			work_data[cont++] = (ixx * grid.ngrid + iyy) * grid.ngrid + izz;
            #endif
         } // izz
        } // iyy
      } // ixx         

      dis = 0.0;
      for(k=Seg[i].start+1;k<Seg[i].start+Seg[i].size;k++)
      {

        delta[0] = Gr[k].Pos[0] - Gr[k-1].Pos[0];
        delta[1] = Gr[k].Pos[1] - Gr[k-1].Pos[1];
        delta[2] = Gr[k].Pos[2] - Gr[k-1].Pos[2];
  
        #ifdef PERIODIC
        if(delta[0]> cp.lbox) delta[0] -= cp.lbox;
        if(delta[0]< 0.0)     delta[0] += cp.lbox;
        if(delta[1]> cp.lbox) delta[1] -= cp.lbox;
        if(delta[1]< 0.0)     delta[1] += cp.lbox;
        if(delta[2]> cp.lbox) delta[2] -= cp.lbox;
        if(delta[2]< 0.0)     delta[2] += cp.lbox;
        #endif
        
        dis += sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);

        ///////////////////////////////////////////////

        Pos_cent[0] = Gr[k].Pos[0];
        Pos_cent[1] = Gr[k].Pos[1];
        Pos_cent[2] = Gr[k].Pos[2];
  
        #ifdef PERIODIC
        if(Pos_cent[0]> cp.lbox) Pos_cent[0] -= cp.lbox;
        if(Pos_cent[0]< 0.0)     Pos_cent[0] += cp.lbox;
        if(Pos_cent[1]> cp.lbox) Pos_cent[1] -= cp.lbox;
        if(Pos_cent[1]< 0.0)     Pos_cent[1] += cp.lbox;
        if(Pos_cent[2]> cp.lbox) Pos_cent[2] -= cp.lbox;
        if(Pos_cent[2]< 0.0)     Pos_cent[2] += cp.lbox;
        #endif
  
      	ixc  = (long)((Pos_cent[0])*(type_real)grid.ngrid*(1.f/cp.lbox));
      	iyc  = (long)((Pos_cent[1])*(type_real)grid.ngrid*(1.f/cp.lbox));
      	izc  = (long)((Pos_cent[2])*(type_real)grid.ngrid*(1.f/cp.lbox));
  
        #ifndef PERIODIC
          for(ixx = ((ixc-1<0) ? 0 : ixc-1); ixx <= ((ixc+1 >= grid.ngrid) ? grid.ngrid-1 : ixc+1); ixx++)
        #else
          for(ixx = ixc-1; ixx <= ixc+1; ixx++)
        #endif
        {
          #ifndef PERIODIC
            for(iyy = ((iyc-1<0) ? 0 : iyc-1); iyy <= ((iyc+1 >= grid.ngrid) ? grid.ngrid-1 : iyc+1); iyy++)
          #else
            for(iyy = iyc-1 ; iyy <= iyc+1 ; iyy++)
          #endif
          {
            #ifndef PERIODIC
              for(izz = ((izc-1<0) ? 0 : izc-1); izz <= ((izc+1 >= grid.ngrid) ? grid.ngrid-1 : izc+1); izz++)
            #else
              for(izz = izc-1 ; izz <= izc+1 ; izz++)
            #endif
            {
  
            	#ifdef PERIODIC
              	work_data[cont++] = ( ( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) )   * grid.ngrid +\
  	    			                        ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ) ) * grid.ngrid +\
                                      ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) );
              #else
          			work_data[cont++] = (ixx * grid.ngrid + iyy) * grid.ngrid + izz;
              #endif
           } // izz
          } // iyy
        } // ixx         
      } // cierra el for de k
  
      qsort(work_data,cont,sizeof(long),&compare_ids);
  
      ixx = 0;   
      for(iyy=0; iyy<cont-1; iyy++) 
        if(work_data[iyy] != work_data[iyy+1]) 
          work_data[ixx++] = work_data[iyy]; 
  
      work_data[ixx++] = work_data[cont-1];   
      cont = ixx;
 
      //dis /= (type_real)Seg[i].size-1;   // LONGITUD PROMEDIO DE SEGMENTOS
      //dis *= M_PI*(rcil2[1] - rcil2[0]); // VOLUMEN DEL PRIMER CILINDRO
      //dis = pow(r_aux,3.0)/dis;          // VOLUMEN DE LA CELDA    
      //izz = NRANDOM*(long)cont*(long)dis;  // CANTIDAD DE RANDOM 
      //naux = izz>(NRANDOM*(long)cp.npart) ? \
      //(NRANDOM*(long)cp.npart) : izz;   // CANTIDAD DE RANDOM si se pasa de un valor critico lo cambio

      //#pragma omp parallel for num_threads(NTHREADS) \
      //schedule(dynamic) default(none) \
      //private(k,j,id_min,rsep,dis,dis_min,naux,\
      //vaux,delta,Pos_cent,iyy,ixx,izz,ixc,iyc,izc,Tid) \
      //shared(cp,ncil,i,cont,RAUX,RINIT,Pos_rand,r_aux,rsep_min,\
      //Seg,Gr,P,grid,dens_data,dens_rand,work_data,stdout) \
      //reduction(+:mu,mass)
      for(iyy=0; iyy<cont; iyy++) 
      {
        Tid = omp_get_thread_num();
        ixx = grid.llirst[work_data[iyy]];
  
        while(ixx != grid.nobj)
        {
           id_min  = -1;
           dis_min = rsep_min = 1e26;
  
           Pos_cent[0] = P[ixx].Pos[0];
           Pos_cent[1] = P[ixx].Pos[1];
           Pos_cent[2] = P[ixx].Pos[2];
  
           // EL ULTIMO NODO
           delta[0] = Pos_cent[0] - Gr[Seg[i].start+Seg[i].size-1].Pos[0];
           delta[1] = Pos_cent[1] - Gr[Seg[i].start+Seg[i].size-1].Pos[1];
           delta[2] = Pos_cent[2] - Gr[Seg[i].start+Seg[i].size-1].Pos[2];
  
           #ifdef PERIODIC
           if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
           if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
           if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
           if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
           if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
           if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
           #endif
  
           dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
           if(dis < Seg[i].Rvir_2[1]) // uso el cuadrado  
           {
             ixx = grid.ll[ixx];  
             continue;
           }
           
           // EL PRIMER NODO
           delta[0] = Pos_cent[0] - Gr[Seg[i].start].Pos[0];
           delta[1] = Pos_cent[1] - Gr[Seg[i].start].Pos[1];
           delta[2] = Pos_cent[2] - Gr[Seg[i].start].Pos[2];
  
           #ifdef PERIODIC
           if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
           if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
           if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
           if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
           if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
           if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
           #endif
  
           dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
           if(dis < Seg[i].Rvir_2[0]) // uso el cuadrado  
           {
             ixx = grid.ll[ixx];  
             continue;
           }

           #ifdef ESFERA
           type_int in_sph_ext = 0;
           #endif
 
           vaux[0] = Gr[Seg[i].start+1].Pos[0]-Gr[Seg[i].start].Pos[0];
           vaux[1] = Gr[Seg[i].start+1].Pos[1]-Gr[Seg[i].start].Pos[1];
           vaux[2] = Gr[Seg[i].start+1].Pos[2]-Gr[Seg[i].start].Pos[2];
  
           #ifdef PERIODIC
           if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
           if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
           if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
           if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
           if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
           if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
           #endif
           rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
           vaux[0] /= rsep;
           vaux[1] /= rsep;
           vaux[2] /= rsep;
  
           // Reutilizo DELTA de lo calculado
           dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];

           // SI ENTRA EN MI CILINDRO
           if(dis>0.0f && dis<rsep) 
           {
             dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;
  
             if(dis>0.0)
             {
               dis = sqrt(dis);
  
               #ifdef BIN_LOG
               dis = log10(dis);
               #endif
  
               if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
               { 
             	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                 if(j < ncil && j >= 0)
                 {
                   id_min = j; // queda igual sumo en el k = 0
                   dis_min = dis;
                   rsep_min = rsep;
                 }
               }
             }
           }

           #ifdef ESFERA
           if(dis>rsep) in_sph_ext = 1;
           #endif

           for(k=Seg[i].start+2;k<Seg[i].start+Seg[i].size;k++)
           {  
             vaux[0] = Gr[k].Pos[0]-Gr[k-1].Pos[0];
             vaux[1] = Gr[k].Pos[1]-Gr[k-1].Pos[1];
             vaux[2] = Gr[k].Pos[2]-Gr[k-1].Pos[2];
  
             #ifdef PERIODIC
             if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
             if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
             if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
             if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
             if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
             if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
             #endif
             rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
             vaux[0] /= rsep;
             vaux[1] /= rsep;
             vaux[2] /= rsep;
  
             delta[0] = Pos_cent[0] - Gr[k-1].Pos[0];
             delta[1] = Pos_cent[1] - Gr[k-1].Pos[1];
             delta[2] = Pos_cent[2] - Gr[k-1].Pos[2];
  
             #ifdef PERIODIC
             if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
             if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
             if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
             if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
             if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
             if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
             #endif
           
             dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];
  
             // SI ENTRA EN MI CILINDRO
             if(dis>0.0f && dis<rsep) 
             {
               dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;
               
               if(dis>0.0)
               {
                 dis = sqrt(dis);
  
                 #ifdef BIN_LOG
                 dis = log10(dis);
                 #endif
  
                 if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                 { 
               	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                   if(j < ncil && j >= 0)
                   {
                     id_min = (k-1-Seg[i].start)*ncil+j;
                     dis_min = dis;                     
                     rsep_min = rsep;
                   }
                 }
               }
  
            #ifdef ESFERA

             in_sph_ext = 0;

             }else if(dis<0.0f){

               if(in_sph_ext==1) 
               {
                 dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
  
                 if(dis>0.0)
                 {
                   dis = sqrt(dis);
  
                   #ifdef BIN_LOG
                   dis = log10(dis);
                   #endif
  
                   if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                   { 
               	    j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                     if(j < ncil && j >= 0)
                     {
                       id_min = (k-1-Seg[i].start)*ncil+j;
                       dis_min = dis;
                       rsep_min = rsep;
                     }
                   }
                 }                           
               }

               in_sph_ext = 0;

             }else{

               in_sph_ext = 1;

             #endif
             }

           } // cierra el for de k
  
           if(id_min != -1)
           {
             k = id_min/ncil;
             j = id_min - ncil*k;   
             j = ncil*Tid + j;
             dens_data[k][j] += 1;

             if(dis_min <= RSPH)
             {
               mu   += 1.0/rsep_min;
               mass += 1;
             }
           }
           
           ixx = grid.ll[ixx];  
  
        } //fin lazo particulas del grid
        
        decode(work_data[iyy],grid.ngrid,&ixc,&iyc,&izc); 

        for(ixx = (iyy == 0      ? 0    : iyy    *DIV_CEIL(naux,(long)cont)); \
            ixx < (iyy == cont-1 ? naux : (iyy+1)*DIV_CEIL(naux,(long)cont)); \
            ixx++)
        {
          id_min  = -1;
          dis_min = 1e26;

          Pos_cent[0] = (Pos_rand[3*ixx]   + (type_real)ixc)*r_aux;
          Pos_cent[1] = (Pos_rand[3*ixx+1] + (type_real)iyc)*r_aux;
          Pos_cent[2] = (Pos_rand[3*ixx+2] + (type_real)izc)*r_aux;  
 
          #ifdef PERIODIC
          if(Pos_cent[0]> cp.lbox) Pos_cent[0] -= cp.lbox;
          if(Pos_cent[0]< 0.0)     Pos_cent[0] += cp.lbox;
          if(Pos_cent[1]> cp.lbox) Pos_cent[1] -= cp.lbox;
          if(Pos_cent[1]< 0.0)     Pos_cent[1] += cp.lbox;
          if(Pos_cent[2]> cp.lbox) Pos_cent[2] -= cp.lbox;
          if(Pos_cent[2]< 0.0)     Pos_cent[2] += cp.lbox;
          #endif
  
          // EL ULTIMO NODO
          delta[0] = Pos_cent[0] - Gr[Seg[i].start+Seg[i].size-1].Pos[0];
          delta[1] = Pos_cent[1] - Gr[Seg[i].start+Seg[i].size-1].Pos[1];
          delta[2] = Pos_cent[2] - Gr[Seg[i].start+Seg[i].size-1].Pos[2];
  
          #ifdef PERIODIC
          if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
          if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
          if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
          if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
          if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
          if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
          #endif
  
          dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
          if(dis < Seg[i].Rvir_2[1]) // uso el cuadrado  
            continue;
          
          // EL PRIMER NODO
          delta[0] = Pos_cent[0] - Gr[Seg[i].start].Pos[0];
          delta[1] = Pos_cent[1] - Gr[Seg[i].start].Pos[1];
          delta[2] = Pos_cent[2] - Gr[Seg[i].start].Pos[2];
  
          #ifdef PERIODIC
          if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
          if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
          if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
          if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
          if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
          if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
          #endif
  
          dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
  
          if(dis < Seg[i].Rvir_2[0]) // uso el cuadrado  
            continue;
  
          vaux[0] = Gr[Seg[i].start+1].Pos[0]-Gr[Seg[i].start].Pos[0];
          vaux[1] = Gr[Seg[i].start+1].Pos[1]-Gr[Seg[i].start].Pos[1];
          vaux[2] = Gr[Seg[i].start+1].Pos[2]-Gr[Seg[i].start].Pos[2];
  
          #ifdef PERIODIC
          if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
          if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
          if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
          if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
          if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
          if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
          #endif
          rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
          vaux[0] /= rsep;
          vaux[1] /= rsep;
          vaux[2] /= rsep;

          #ifdef ESFERA
          type_int in_sph_ext = 0;
          #endif

          // Reutilizo DELTA de lo calculado
          dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];
  
          // SI ENTRA EN MI CILINDRO
          if(dis>0.0f && dis<rsep) 
          {
            dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;
  
            if(dis>0.0)
            {
              dis = sqrt(dis);
  
              #ifdef BIN_LOG
              dis = log10(dis);
              #endif
  
              if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
              { 
            	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
  
                if(j < ncil && j >= 0)
                {
                  id_min = j; // queda igual sumo en el k = 0
                  dis_min = dis;
                }
              }
            }
          }

          #ifdef ESFERA
          if(dis>rsep) in_sph_ext = 1;
          #endif

          for(k=Seg[i].start+2;k<Seg[i].start+Seg[i].size;k++)
          {  
            vaux[0] = Gr[k].Pos[0]-Gr[k-1].Pos[0];
            vaux[1] = Gr[k].Pos[1]-Gr[k-1].Pos[1];
            vaux[2] = Gr[k].Pos[2]-Gr[k-1].Pos[2];
  
            #ifdef PERIODIC
            if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
            if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
            if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
            if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
            if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
            if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
            #endif
            rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
            vaux[0] /= rsep;
            vaux[1] /= rsep;
            vaux[2] /= rsep;
  
            delta[0] = Pos_cent[0] - Gr[k-1].Pos[0];
            delta[1] = Pos_cent[1] - Gr[k-1].Pos[1];
            delta[2] = Pos_cent[2] - Gr[k-1].Pos[2];
  
            #ifdef PERIODIC
            if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
            if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
            if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
            if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
            if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
            if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
            #endif
          
            dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];
  
            // SI ENTRA EN MI CILINDRO
            if(dis>0.0f && dis<rsep) 
            {
              dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;

              if(dis>0.0)
              {
                dis = sqrt(dis);
  
                #ifdef BIN_LOG
                dis = log10(dis);
                #endif
  
                if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                { 
              	  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
  
                  if(j < ncil && j >= 0)
                  {
                    id_min = (k-1-Seg[i].start)*ncil+j;
                    dis_min = dis;
                  }
                }
              }
            
            #ifdef ESFERA
             in_sph_ext = 0;

            }else if(dis<0.0f){

              if(in_sph_ext==1) 
              {
                dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
  
                if(dis>0.0)
                {
                  dis = sqrt(dis);
  
                  #ifdef BIN_LOG
                  dis = log10(dis);
                  #endif
  
                  if((dis > RINIT) && (dis < RAUX) && (dis<dis_min))
                  { 
              	    j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
      
                    if(j < ncil && j >= 0)
                    {
                      id_min = (k-1-Seg[i].start)*ncil+j;
                      dis_min = dis;
                    }
                  }
                }                           
              }                           

              in_sph_ext = 0;

             }else{

               in_sph_ext = 1;

              #endif
            }
  
          } // cierra el for de k
  
          if(id_min != -1)
          {
            k = id_min/ncil;
            j = id_min - ncil*k;   
            j = ncil*Tid + j;
            dens_rand[k][j] += 1;
          }
  
        } // cierra el loop sobre las celdas

      } // cierra la cantidad de random      
      
      rho   = (type_real)cont*pow(r_aux,3.0); // vol tot
      rho   = (type_real)naux/rho;            // dens tot
      //rho  /= (type_real)(Seg[i].size-1);     // normalizacion por segmento
  
      cont = 0;
      for(k=0;k<Seg[i].size-1;k++)
      {
        for(j=ncil;j<NTHREADS*ncil;j++)
        {
          dens_data[k][j%ncil] += dens_data[k][j];
          dens_rand[k][j%ncil] += dens_rand[k][j];
        }

        id_min = 1;
        for(j=0;j<ncil;j++)
        {
          id_min *= (dens_rand[k][j]>0 ? 1 : 0);
          tmp_matrix_mean_cil[j]      += (type_real)dens_data[k][j];
          tmp_matrix_mean_cil[ncil+j] += (type_real)dens_rand[k][j];
        }
  
        if(id_min == 1)
        {
          for(j=0;j<ncil;j++)
            tmp_matrix_mean_profile[j] += rho*((type_real)dens_data[k][j]/(type_real)dens_rand[k][j]);
          cont++;
        }
  
        free(dens_rand[k]);
        free(dens_data[k]);
      }
  
      pf = fopen(filename,"w");
  
      #ifdef BINARIO
        fwrite(&ncil,sizeof(type_int),1,pf);        
      #endif
  
      for(j=0;j<ncil;j++)
      {
        dis = 0.5*(sqrt(rcil2[j])+sqrt(rcil2[j+1]))/1000.0; // En Mpc
        tmp_matrix_mean_profile[j] = tmp_matrix_mean_profile[j]*pow(1000.0,3);           // En Mpc
      	tmp_matrix_mean_profile[j] /= (type_real)(cont);     // normalizacion por segmento
  
        tmp_matrix_mean_cil[j]  = rho*(tmp_matrix_mean_cil[j]/tmp_matrix_mean_cil[ncil+j]);
        tmp_matrix_mean_cil[j] *= pow(1000.0,3);           // En Mpc
        // tmp_matrix_mean_cil[j] /= (type_real)(Seg[i].size-1);
  
        #ifdef BINARIO
          fwrite(&dis,sizeof(type_real),1,pf);        
          fwrite(&tmp_matrix_mean_profile[j],sizeof(type_real),1,pf);        
        #else
          fprintf(pf,"%d %d %f %f %f\n",cont,Seg[i].size-1,dis,tmp_matrix_mean_profile[j],tmp_matrix_mean_cil[j]);
        #endif
  
      }
      fclose(pf);
 
      set_name_prop(filename,Seg[i].id);
      pf = fopen(filename,"w");

      mu  *= (cp.Mpart*1000.0);        // En Mass[10^10 Msol / h]/Mpc
      rsep = cp.Mpart*(type_real)mass; // En Mass[10^10 Msol / h]
      dis  = mu/(pow(RSPH/1000.0,2));  // En Mass[10^10 Msol / h]/Mpc^3

      #ifdef BINARIO
        fwrite(&rsep,sizeof(type_real),1,pf);        
        fwrite(&mu,sizeof(type_real),1,pf);        
        fwrite(&dis,sizeof(type_real),1,pf);        
      #else
        fprintf(pf,"%f %f %f\n",rsep,mu,dis);
      #endif
      fclose(pf);

      free(work_data);
      free(dens_rand);
      free(dens_data);
      free(tmp_matrix_mean_cil);
      free(tmp_matrix_mean_profile);
  
    }// cierro el for sobre los filamentos
  
  ////////////////////////////////////////////////////////////////////////////////////////////  
    free(Pos_rand);

    return;
  
  }

#endif

#else

  static void calc_profile(void)
  {
    int i, j, k;
    long ixx, iyy, izz;
    long ixc, iyc, izc, naux;
    type_real dis, rsep;
    type_real vaux[3], delta[3], Pos_cent[3];
    type_real *tmp_matrix_mean_profile;
    type_real *tmp_matrix_mu;
    type_int  *tmp_matrix_mass;
    FILE *pf;
    char filename[200];
#ifdef ROTACION
    type_real rot_matrix[9];
    type_real vrot[3];
#endif
    const type_int NTHREADS = omp_get_max_threads();

    for(i=0;i<cp.nseg;i++)
    {
  
      tmp_matrix_mean_profile = (type_real  *) calloc(NTHREADS*ncil,sizeof(type_real));
      tmp_matrix_mass         = (type_int   *) calloc(NTHREADS,sizeof(type_int));
      tmp_matrix_mu          = (type_real  *) calloc(NTHREADS,sizeof(type_real));
  
      fprintf(stdout,"size %d seg %d/%d fil %d/%d\n",Seg[i].size-1,Seg[i].start-i,cp.ngrup-cp.nseg,i,cp.nseg);
      fflush(stdout);
  
#ifdef ROTACION
      #pragma omp parallel for\
      num_threads(NTHREADS) default(none) \
      private(k,j,dis,vaux,delta,Pos_cent,\
      ixx,iyy,izz,ixc,iyc,izc,naux, \
      rot_matrix,vrot) \
      shared(cp,ncil,i,RAUX,RINIT,rsep,Seg,Gr,P,grid,stdout,\
      tmp_matrix_mean_profile,tmp_matrix_mu,tmp_matrix_mass) 
#else
      #pragma omp parallel for\
      num_threads(NTHREADS) default(none) \
      private(k,j,dis,vaux,delta,Pos_cent,\
      ixx,iyy,izz,ixc,iyc,izc,naux) \
      shared(cp,ncil,i,RAUX,RINIT,rsep,Seg,Gr,P,grid,stdout,\
      tmp_matrix_mean_profile,tmp_matrix_mu,tmp_matrix_mass) 
#endif
      for(k=Seg[i].start+1;k<Seg[i].start+Seg[i].size;k++)
      {
  
        type_int  Tid = omp_get_thread_num(); 
  
        vaux[0] = Gr[k].Pos[0]-Gr[k-1].Pos[0];
        vaux[1] = Gr[k].Pos[1]-Gr[k-1].Pos[1];
        vaux[2] = Gr[k].Pos[2]-Gr[k-1].Pos[2];

        #ifdef PERIODIC
        if(vaux[0]> 0.5*cp.lbox) vaux[0] -= cp.lbox;
        if(vaux[0]<-0.5*cp.lbox) vaux[0] += cp.lbox;
        if(vaux[1]> 0.5*cp.lbox) vaux[1] -= cp.lbox;
        if(vaux[1]<-0.5*cp.lbox) vaux[1] += cp.lbox;
        if(vaux[2]> 0.5*cp.lbox) vaux[2] -= cp.lbox;
        if(vaux[2]<-0.5*cp.lbox) vaux[2] += cp.lbox;
        #endif

        rsep = sqrt(vaux[0]*vaux[0]+vaux[1]*vaux[1]+vaux[2]*vaux[2]);
  
        vaux[0] /= rsep;
        vaux[1] /= rsep;
        vaux[2] /= rsep;
  
        Pos_cent[0] = 0.5*rsep*vaux[0] + Gr[k-1].Pos[0];
        Pos_cent[1] = 0.5*rsep*vaux[1] + Gr[k-1].Pos[1];
        Pos_cent[2] = 0.5*rsep*vaux[2] + Gr[k-1].Pos[2];
  
        #ifdef PERIODIC
        if(Pos_cent[0]> cp.lbox) Pos_cent[0] -= cp.lbox;
        if(Pos_cent[0]< 0.0)     Pos_cent[0] += cp.lbox;
        if(Pos_cent[1]> cp.lbox) Pos_cent[1] -= cp.lbox;
        if(Pos_cent[1]< 0.0)     Pos_cent[1] += cp.lbox;
        if(Pos_cent[2]> cp.lbox) Pos_cent[2] -= cp.lbox;
        if(Pos_cent[2]< 0.0)     Pos_cent[2] += cp.lbox;
        #endif
  
        ixc  = (long)((Pos_cent[0])*(type_real)grid.ngrid*(1.f/cp.lbox));
        iyc  = (long)((Pos_cent[1])*(type_real)grid.ngrid*(1.f/cp.lbox));
        izc  = (long)((Pos_cent[2])*(type_real)grid.ngrid*(1.f/cp.lbox));
  
    	  #ifndef PERIODIC
    	    for(ixx = ((ixc-1<0) ? 0 : ixc-1); ixx <= ((ixc+1 >= grid.ngrid) ? grid.ngrid-1 : ixc+1); ixx++)
    	  #else
    	    for(ixx = ixc-1; ixx <= ixc+1; ixx++)
    	  #endif
    	  {
    	    #ifndef PERIODIC
    	      for(iyy = ((iyc-1<0) ? 0 : iyc-1); iyy <= ((iyc+1 >= grid.ngrid) ? grid.ngrid-1 : iyc+1); iyy++)
    	    #else
    	      for(iyy = iyc-1 ; iyy <= iyc+1 ; iyy++)
    	    #endif
    	    {
    	      #ifndef PERIODIC
    	        for(izz = ((izc-1<0) ? 0 : izc-1); izz <= ((izc+1 >= grid.ngrid) ? grid.ngrid-1 : izc+1); izz++)
    	      #else
    	        for(izz = izc-1 ; izz <= izc+1 ; izz++)
    	      #endif
    	      {
  
    	      	#ifdef PERIODIC
    	        	naux = ( ( (ixx >= (long)grid.ngrid) ? ixx-(long)grid.ngrid : ( (ixx<0) ? ixx + (long)grid.ngrid : ixx ) )   * grid.ngrid +\
  		  			            ( (iyy >= (long)grid.ngrid) ? iyy-(long)grid.ngrid : ( (iyy<0) ? iyy + (long)grid.ngrid : iyy ) ) ) * grid.ngrid +\
    	                    ( (izz >= (long)grid.ngrid) ? izz-(long)grid.ngrid : ( (izz<0) ? izz + (long)grid.ngrid : izz ) );
    	        #else
    	    			naux = (ixx * grid.ngrid + iyy) * grid.ngrid + izz;
    	        #endif
  
              naux = grid.llirst[naux];
  
              while(naux != grid.nobj)
              {

                // EL ULTIMO NODO
                delta[0] = P[naux].Pos[0] - Gr[Seg[i].start+Seg[i].size-1].Pos[0];
                delta[1] = P[naux].Pos[1] - Gr[Seg[i].start+Seg[i].size-1].Pos[1];
                delta[2] = P[naux].Pos[2] - Gr[Seg[i].start+Seg[i].size-1].Pos[2];

                #ifdef PERIODIC
                if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
                if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
                if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
                if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
                if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
                if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
                #endif

                dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];

                if(dis < Seg[i].Rvir_2[1]) // uso el cuadrado  
                {
                  naux = grid.ll[naux];  
                  continue;
                }
                
                // EL PRIMER NODO
                delta[0] = P[naux].Pos[0] - Gr[Seg[i].start].Pos[0];
                delta[1] = P[naux].Pos[1] - Gr[Seg[i].start].Pos[1];
                delta[2] = P[naux].Pos[2] - Gr[Seg[i].start].Pos[2];

                #ifdef PERIODIC
                if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
                if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
                if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
                if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
                if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
                if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
                #endif

                dis = delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];

                if(dis < Seg[i].Rvir_2[0]) // uso el cuadrado  
                {
                  naux = grid.ll[naux];  
                  continue;
                }

                delta[0] = P[naux].Pos[0] - Gr[k-1].Pos[0];
                delta[1] = P[naux].Pos[1] - Gr[k-1].Pos[1];
                delta[2] = P[naux].Pos[2] - Gr[k-1].Pos[2];
  
                #ifdef PERIODIC
                if(delta[0]> 0.5*cp.lbox) delta[0] -= cp.lbox;
                if(delta[0]<-0.5*cp.lbox) delta[0] += cp.lbox;
                if(delta[1]> 0.5*cp.lbox) delta[1] -= cp.lbox;
                if(delta[1]<-0.5*cp.lbox) delta[1] += cp.lbox;
                if(delta[2]> 0.5*cp.lbox) delta[2] -= cp.lbox;
                if(delta[2]<-0.5*cp.lbox) delta[2] += cp.lbox;
                #endif
  
              #ifdef ROTACION 

                dis  = sqrt(vaux[1]*vaux[1]+vaux[0]*vaux[0]);
                vrot[0] = -vaux[0]/dis;  //sin_phi   = -vaux[0]/dis;
                vrot[1] =  vaux[1]/dis;  //cos_phi   =  vaux[1]/dis;
                vrot[2] =  vaux[2];      //cos_theta =  vaux[2]/rsep; // no lo dividio por el modulo
                dis     = -dis;          //sin_theta = -dis/rsep;     //

                rot_matrix[0] =  vrot[1]         ; rot_matrix[1] =  vrot[0]         ; rot_matrix[2] = 0.0;
                rot_matrix[3] = -vrot[0]*vrot[2] ; rot_matrix[4] =  vrot[2]*vrot[1] ; rot_matrix[5] = dis;
                rot_matrix[6] =  dis*vrot[0]     ; rot_matrix[7] = -dis*vrot[1]     ; rot_matrix[8] = vrot[2];

            	  for(j=0; j<3; j++)
                  vrot[j] = rot_matrix[3*j]*delta[0] + rot_matrix[3*j+1]*delta[1] + rot_matrix[3*j+2]*delta[2];

                // SI ENTRA EN MI CILINDRO
                if(vrot[2]>0.0 && vrot[2]<rsep)
                {
                  dis  = sqrt(vrot[1]*vrot[1]+vrot[0]*vrot[0]);

            #else

                dis = vaux[0]*delta[0] + vaux[1]*delta[1] + vaux[2]*delta[2];

                // SI ENTRA EN MI CILINDRO
                if(dis>0.0f && dis<rsep) 
                {
                  dis = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] - dis*dis;

                  if(dis>0.0)
                  {

                    dis = sqrt(dis);
            #endif
              
                  if(dis <= RSPH)
                  {
                    tmp_matrix_mu[Tid]   += 1.0/rsep;
                    tmp_matrix_mass[Tid] += 1;
                  }

                  #ifdef BIN_LOG
                  dis = log10(dis);
                  #endif
 
                  if((dis <= RINIT) || (dis >= RAUX))
                  {
                    naux = grid.ll[naux];                          
                    continue;
                  }
  
                  j = (int)((dis-RINIT)/(RAUX-RINIT)*ncil);
  
                  if(j >= ncil || j <= -1)
                  {
                    naux = grid.ll[naux];                          
                    continue;
                  }
  
                  j = ncil*Tid+j;
                  tmp_matrix_mean_profile[j] += 1.0/rsep;
              #ifndef ROTACION
                }
              #endif

                }
  
                naux = grid.ll[naux];
  
              } //fin lazo particulas del grid
             } // izz
            } // iyy
          } // ixx         
      }// cierro el for sobre los segmentos
  
      set_name(filename,Seg[i].id);

      pf = fopen(filename,"w");

      #ifdef BINARIO
        fwrite(&ncil,sizeof(type_int),1,pf);        
      #endif

      for(j=ncil;j<NTHREADS*ncil;j++)
        tmp_matrix_mean_profile[j%ncil] += tmp_matrix_mean_profile[j];
  
      for(j=0;j<ncil;j++)
      {
        tmp_matrix_mean_profile[j] /= M_PI*(rcil2[j+1]-rcil2[j]);
        tmp_matrix_mean_profile[j] /= (type_real)(Seg[i].size-1);
  
        dis = 0.5*(sqrt(rcil2[j])+sqrt(rcil2[j+1]))/1000.0; // En Mpc
        tmp_matrix_mean_profile[j] = tmp_matrix_mean_profile[j]*pow(1000.0,3);           // En Mpc
  
        #ifdef BINARIO
            fwrite(&dis,sizeof(type_real),1,pf);        
            fwrite(&tmp_matrix_mean_profile[j],sizeof(type_real),1,pf);        
        #else
            fprintf(pf,"%f %f\n",dis,tmp_matrix_mean_profile[j]);
        #endif
      }
      fclose(pf);

      set_name_prop(filename,Seg[i].id);
      pf = fopen(filename,"w");

      for(j=1;j<NTHREADS;j++)
      {
        tmp_matrix_mu[0]   += tmp_matrix_mu[j];
        tmp_matrix_mass[0] += tmp_matrix_mass[j];
      }
      
      tmp_matrix_mu[0] *= (cp.Mpart*1000.0);         // En Mass[10^10 Msol / h]/Mpc
      rsep = cp.Mpart*(type_real)tmp_matrix_mass[0]; // En Mass[10^10 Msol / h]
      dis  = tmp_matrix_mu[0]/(pow(RSPH/1000.0,2));  // En Mass[10^10 Msol / h]/Mpc^3

      #ifdef BINARIO
        fwrite(&rsep,sizeof(type_real),1,pf);        
        fwrite(&tmp_matrix_mu[0],sizeof(type_real),1,pf);        
        fwrite(&dis,sizeof(type_real),1,pf);        
      #else
        fprintf(pf,"%f %f %f\n",rsep,tmp_matrix_mu[0],dis);
      #endif
      fclose(pf);

      free(tmp_matrix_mu);
      free(tmp_matrix_mass);
      free(tmp_matrix_mean_profile);
    }// cierro el for sobre los filamentos
  
    ////////////////////////////////////////////////////////////////////////////////////////////  
  
    return;
  }

#endif

extern void propiedades(type_real *fof)
{

  r_aux = sqrt(2.0f)*RAUX; // RADIO DE BUSQUEDA EN Kpc

#ifdef BIN_LOG
    RAUX  = log10(RAUX);
    RINIT = log10(RINIT);  
#endif 
 
  rcil2  = (type_real *) malloc((ncil+1)*sizeof(type_real));
  #ifdef BIN_LOG
    logspace(rcil2,RAUX,RINIT,(ncil+1));
  #else
    linspace(rcil2,RAUX,RINIT,(ncil+1));
  #endif
  
  grid.nobj = cp.npart;
  grid.ngrid = (long)(cp.lbox/r_aux);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
    r_aux = cp.lbox/(type_real)NGRIDMAX;
  }else{
    fprintf(stdout,"Using NGRID = %lu\n",grid.ngrid);
  }

  fflush(stdout);

  grid_init();
  grid_build();

  BLUE("******************************\n");
  GREEN("******************************\n");
  fflush(stdout);

#ifdef SIN_REPETICION
  calc_profile_sin_repeticion();
#else
  calc_profile();
#endif

  GREEN("******************************\n");
  BLUE("******************************\n");
  fflush(stdout);

  grid_free();

  free(rcil2);

  return;
}
