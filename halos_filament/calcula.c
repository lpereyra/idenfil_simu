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

void propiedades(int NNN, type_real *fof)
{
  int totfin = 0;
  type_int i, k, id;
  type_real RMAX, RMIN, LMAX;
  type_real rsep;                   // En Kpc
  type_real rlong;                  // En Kpc
  const type_real rcil = 10.0f*1000.0f; // RADIO DE BUSQUEDA EN Kpc
  char filename[200], name[200];
  FILE *pf;
  #ifdef MPC 
  type_real POSFACTOR = 1000.;
  #endif
    
  #ifdef MPC
  GREEN("********** IMPORTANTE ***********\n");
  cp.lbox *= POSFACTOR;
  sprintf(filename,"Reescala lbox %g Mpc to %g Kpc\n",cp.lbox/POSFACTOR,cp.lbox);GREEN(filename);
  GREEN("**********************************\n");
  #endif

  stadistic(cp.nseg,&RMAX,&RMIN,&LMAX); 

  #ifdef FIXED_SEPARATION
    ncil = (int)(RLEN/RSEP);
    RLEN *= 1000.f;
    RSEP *= 1000.f;
    assert(RLEN<RMAX);
    sprintf(message,"Num cylindres:          %d\n",ncil);BLUE(message);
  #endif

  fprintf(stdout,"RMAX %f Mpc\n",RMAX/1000.);
  fprintf(stdout,"RMIN %f Mpc\n",RMIN/1000.);
  fprintf(stdout,"LMAX %f Mpc\n",LMAX/1000.);

  grid.nobj = cp.ngrup;
  grid.ngrid = (int)(cp.lbox/2000.);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }else{
    fprintf(stdout,"Using NGRID = %lu\n",grid.ngrid);
  }

  grid_init();
  grid_build();

  BLUE("******************************\n");
  GREEN("******************************\n");

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) private(i,rsep,rlong) \
  shared(Seg,cp,nbins)
  for(i=0;i<cp.nseg;i++)
  {
    rsep    = Seg[i].len/(type_real)(nbins-1);
    rlong   = 2.0*rsep; 

    calcular_pert(i,rsep,rlong,0.5f*rlong,rcil*rcil);
  } //FINALIZA EL PARALELO
  BLUE("******************************\n");

  set_name(name,NNN,"halos_por_filament");
  sprintf(filename,"%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pf = fopen(filename,"w");

  fwrite(&cp.nseg,sizeof(int),1,pf);        

  for(i=0;i<cp.nseg;i++)
  {
    fwrite(&i,sizeof(int),1,pf);                   // Num Filamento
    fwrite(&Seg[i].flag,sizeof(int),1,pf);         // escribe la bandera
    fwrite(&Seg[i].len,sizeof(type_real),1,pf);    // escribe la longitud del segmento    

    rsep = Gr[Seg[i].list[0]].NumPart*cp.Mpart;
    fwrite(&rsep,sizeof(type_real),1,pf);             // escribe la masa del primer nodo
    rsep = Gr[Seg[i].list[Seg[i].size-1]].NumPart*cp.Mpart;
    fwrite(&rsep,sizeof(type_real),1,pf);             // escribe la masa del ultimo nodo

    /////////////////////////////////////////////////////////////////////////////////
    k = Seg[i].nlista_halos;
    qsort(Seg[i].lista_halos, k, sizeof(struct data_list), cmp);
    fwrite(&k,sizeof(type_int),1,pf);                           // Num de halos
    fwrite(&Seg[i].lista_halos,sizeof(struct data_list),k,pf);  // Struct halos
    /////////////////////////////////////////////////////////////////////////////////
    totfin++;
  }

  fclose(pf);

  GREEN("******************************\n");
  BLUE("******************************\n");
  fprintf(stdout,"Escribe %d Filamentos\n",totfin);  
  BLUE("******************************\n");

  GREEN("******************************\n");
  BLUE("Free...\n");
  grid_free();
  GREEN("******************************\n");

  totfin = 0;

  for(i=0;i<cp.ngrup;i++)
    Gr[i].nlista_fil = 0;

  for(i=0;i<cp.nseg;i++)
  {
    for(k=0;k<Seg[i].nlista_halos;k++)
    { 
      struct data_list aux_dat = Seg[i].lista_halos[k];
      assert(aux_dat.id<cp.ngrup);
      Gr[aux_dat.id].nlista_fil += 1;
    }  
  }
  
  for(i=0;i<cp.ngrup;i++)
  {
    Gr[i].lista_fil  = (struct data_list *) malloc(Gr[i].nlista_fil*sizeof(struct data_list));
    Gr[i].nlista_fil = 0;
  }
  
  BLUE("******************************\n");
  
  for(i=0;i<cp.nseg;i++)
  {
    for(k=0;k<Seg[i].nlista_halos;k++)
    {
      id = Seg[i].lista_halos[k].id;
      Gr[id].lista_fil[Gr[id].nlista_fil]    = Seg[i].lista_halos[k];
      Gr[id].lista_fil[Gr[id].nlista_fil].id = i;
      Gr[id].nlista_fil++;
    }
  }

  BLUE("******************************\n");

  set_name(name,NNN,"filament_por_halos");
  sprintf(filename,"%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pf = fopen(filename,"w");

  fwrite(&cp.ngrup,sizeof(type_int),1,pf);        
  
  for(i=0;i<cp.ngrup;i++)
  {
    fwrite(&Gr[i].Pos[0],sizeof(type_real),1,pf);
    fwrite(&Gr[i].Pos[1],sizeof(type_real),1,pf);
    fwrite(&Gr[i].Pos[2],sizeof(type_real),1,pf);
    fwrite(&Gr[i].NumPart,sizeof(type_int),1,pf);

    k = Gr[i].nlista_fil;
    qsort(Gr[i].lista_fil, k, sizeof(struct data_list), cmp);

    fwrite(&k,sizeof(type_int),1,pf);                        // Num de fil

    if(k!=0) // por si hay halos sin filamentos
      fwrite(&Gr[i].lista_fil,sizeof(struct data_list),k,pf);  // Struct fil

    totfin++;
  }

  BLUE("******************************\n");
  fprintf(stdout,"Escribe %d Halos\n",totfin);  
  BLUE("******************************\n");

  fclose(pf);

  return;
}

void set_name(char * name, const int NNN, char * prefix)
{

  //#ifdef FIXED_SEPARATION
  //  sprintf(name,"%.4d_fixed",NNN);
  //#else
  //  sprintf(name,"%.4d_varia",NNN);
  //#endif

  sprintf(name,"%.4d",NNN);

  #ifdef MCRITIC
    sprintf(name,"%s_%s_cut",name,prefix);
  #else
    sprintf(name,"%s_%s",name,prefix);
  #endif

  //#ifdef BIN_LOG
  //  sprintf(name,"%s_LOG",name);
  //#else
  //  sprintf(name,"%s_LIN",name);
  //#endif

}

int point_inside(type_real dot, type_real rlong_2)
{
  if(fabs(dot)>rlong_2){

      return 0;

  }else{

    return 1;

  }
}

int cmp(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

void stadistic(int n, type_real *MAX, type_real *MIN, type_real *LMAX) 
{
  int i, j;
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
void logspace(type_real *rcil2, type_real max, type_real min, int bins) 
{

  int i;
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

void linspace(type_real *rcil2, type_real max, type_real min, int bins) 
{

  int i;
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

void calcular_pert(const int i, const type_real rsep, \
         const type_real rlong, const type_real rlong_2, type_real const rcil2)
{

  int k,l,idim;
  int bin,lbin;
  type_real vdir[3];
  type_real Pos_cent[3];
  type_real mod_v[3];
  type_real r,rbin,racum;  
  struct list *aux_lista_halos;

  l = 0;
  racum   = 0.0f;
  aux_lista_halos = NULL;

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
    vdir[idim] *= (1./r);

	////////////////////////////////////////////////////////////////

  calc_cil(Gr[Seg[i].list[0]].Pos,vdir,rcil2,rlong_2,&aux_lista_halos);             
  l++; // sumo

	////////////////////////////////////////////////////////////////

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

 	  ////////////////////////////////////////////////////////////////
    calc_cil(Pos_cent,vdir,rcil2,rlong_2,&aux_lista_halos);             
    l++; // sumo
  	////////////////////////////////////////////////////////////////

    r = 0.0f;
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

 	    ///////////////////////////////////////////////////////
      calc_cil(Pos_cent,vdir,rcil2,rlong_2,&aux_lista_halos);    
      l++; // sumo
      ///////////////////////////////////////////////////////

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

 	////////////////////////////////////////////////////////////////
  calc_cil(Gr[Seg[i].list[Seg[i].size-1]].Pos,vdir,rcil2,rlong_2,&aux_lista_halos);
  l++; // sumo
 	////////////////////////////////////////////////////////////////
  
  assert(l==nbins);

  /////print_lista(aux_lista_halos); 
  remove_duplicates(aux_lista_halos);
  /////fprintf(stdout,"\nsin duplicados\n");
  /////print_lista(aux_lista_halos);

  count_list(aux_lista_halos, &Seg[i].nlista_halos);

  Seg[i].lista_halos = (struct data_list *) malloc(Seg[i].nlista_halos * sizeof(struct data_list));  

  alocate_list(aux_lista_halos,Seg[i].lista_halos);

  free_list(aux_lista_halos);
  
  return;
}

void calc_cil(const type_real * restrict Pos_cent,  const type_real * restrict versor, const type_real rcil2, const type_real rlong_2, struct list ** lista)
{
	
	int  i;
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox,ilim;
  type_real lbox,fac,lbox2;
  type_real dot,dis;
  type_real Posprima[3];

  lbox  = cp.lbox;
  fac   = (type_real)grid.ngrid/lbox;
  lbox2 = lbox/2.0;

  dis = sqrt(rcil2);
  ilim = dis>rlong_2 ? (long)(dis*fac)+1 : (long)(rlong_2*fac)+1;

  ixc  = (int)(Pos_cent[0]*fac);
  ixci = ixc - ilim;
  ixcf = ixc + ilim;
  iyc  = (int)(Pos_cent[1]*fac);
  iyci = iyc - ilim;
  iycf = iyc + ilim;
  izc  = (int)(Pos_cent[2]*fac);
  izci = izc - ilim;
  izcf = izc + ilim;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= grid.ngrid ) ixcf = grid.ngrid - 1;
  if( iycf >= grid.ngrid ) iycf = grid.ngrid - 1;
  if( izcf >= grid.ngrid ) izcf = grid.ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= (long)grid.ngrid) ix = ix - (long)grid.ngrid;
    if(ix < 0) ix = ix + grid.ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= (long)grid.ngrid) iy = iy - (long)grid.ngrid;
      if(iy < 0) iy = iy + (long)grid.ngrid;
      #endif
  
      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= (long)grid.ngrid) iz = iz - (long)grid.ngrid;
        if(iz < 0) iz = iz + (long)grid.ngrid;
        #endif

        ibox = (ix * (long)grid.ngrid + iy) * (long)grid.ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {
          Posprima[0] = Gr[i].Pos[0] - Pos_cent[0];
          Posprima[1] = Gr[i].Pos[1] - Pos_cent[1];
          Posprima[2] = Gr[i].Pos[2] - Pos_cent[2];

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

            if(dis<rcil2)
            {
              struct data_list dat;
              dat.id = i;
              dat.r  = sqrt(fabs(dis)); 
              // le tomo el valor absoluto por el error en la precision
              // en los nodos que estan muy cerca del eje puede dar negativo
              push_list(lista, dat);
            }// cierra el if

          }

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}
