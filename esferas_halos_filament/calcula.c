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

#include <iterator>
#include <utility>
#include <set>


struct info {
  type_int  id;
  type_int  k;
  type_real r_sph;
};

struct data {

  type_int      nfil;
  struct info * matrix;

};

static omp_lock_t *lock; 

static void realloc_set(std::set<std::pair<type_real,type_int> > &myset)
{
  std::set<std::pair<type_real,type_int> >::iterator it = myset.begin(); 
  std::advance (it,K_NEIGH_SIZE);
  myset.erase(it,myset.end());

  return;
}

static void set_name(char * name, const type_int NNN, const char * prefix)
{

  #ifdef ORIGINAL
    sprintf(name,"%.2d_%.4d",snap.num,NNN);
  #else
    sprintf(name,"%.2d_%.4d_smooth",snap.num,NNN);
  #endif

  #ifdef NEW
    sprintf(name,"%s_new",name);
  #endif

  #ifdef MCRITIC
    sprintf(name,"%s_cut_%.2f",name,m_critica);
  #endif

  sprintf(name,"%s_%s",name,prefix);

  #ifdef EXTEND
    sprintf(name,"%s_extend",name);
  #endif

  sprintf(name,"%s_k_neigh_%.2d",name,K_NEIGH_SIZE);

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);

  return;
}

static void calc_search_fil(const type_real * __restrict__ Pos_cent, const type_real r_2, struct list ** lista)
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
          Posprima[0] = Seg_centr[i].Pos[0] - Pos_cent[0];
          Posprima[1] = Seg_centr[i].Pos[1] - Pos_cent[1];
          Posprima[2] = Seg_centr[i].Pos[2] - Pos_cent[2];

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
            dat.id = Seg_centr[i].save;
            dat.k = Seg_centr[i].id;
            dat.r = sqrt(dis);
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
  const type_real rcil2 = rbus*rbus;
  struct list *aux_lista_fil = NULL;
  struct list *root;
  std::pair<type_real,type_int> result;
  type_int k,idfil;

  #ifdef PARTICLES
  calc_search_fil(P[i].Pos,rcil2,&aux_lista_fil);             
  #else
  calc_search_fil(Gr[i].Pos,rcil2,&aux_lista_fil);             
  #endif

  //remove_duplicates(aux_lista_fil);
  
  while(aux_lista_fil != NULL)
  {
    root = aux_lista_fil;
    
    idfil = root->data.id;
    k = root->data.k;
    result = std::make_pair(root->data.r,i);

    omp_set_lock(&lock[idfil]);

    if(Seg[idfil].k_neigh[k].size()>SET_MAX_SIZE)
      realloc_set(Seg[idfil].k_neigh[k]);

    Seg[idfil].k_neigh[k].insert(result);

    omp_unset_lock(&lock[idfil]);

    aux_lista_fil = aux_lista_fil->next;
    free(root);
  }

  return;
}

extern void propiedades(type_int NNN, type_real *fof, const type_real r_aux)
{
  type_int i,j,k;
  type_real dens, mass, len;
  std::pair<type_real, type_int> rs;
  char filename[200];
  FILE *pfout, *pfpropiedades;


  lock = (omp_lock_t *) malloc(cp.nseg*sizeof(omp_lock_t));

  len = 0.0f;
  for(i=0;i<cp.nseg;i++)
  {
   len += (Seg[i].size*SET_MAX_SIZE)*sizeof(rs);
   omp_init_lock(&(lock[i]));
  }

  fprintf(stdout,"allocating max %.5f gb\n",
      (double)(len/1024.0/1024.0/1024.0));

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  BLUE("******************************\n");
  GREEN("******************************\n");
  fflush(stdout);

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) \
  private(i) shared(cp,stdout)
  #ifdef PARTICLES
    for(i=0;i<cp.npart;i++)
  #else
    for(i=0;i<cp.ngrup;i++)
  #endif
  {
    #ifdef PARTICLES
    if(i%1000000==0) fprintf(stdout,"%u %.4f\n",i,(type_real)i/(type_real)cp.npart);
    #else
    if(i%10000==0) fprintf(stdout,"%u %.4f\n",i,(type_real)i/(type_real)cp.ngrup);
    #endif
    calcular_pert(i,r_aux);
  }

  set_name(filename,NNN,"segmentos");
  pfout = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pfout);

  set_name(filename,NNN,"propiedades");
  pfpropiedades=fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(type_int),1,pfpropiedades);

  for(i=0;i<cp.nseg;i++)
  {   
    fwrite(&Seg[i].flag,sizeof(type_int),1,pfpropiedades);
    fwrite(&Seg[i].size,sizeof(type_int),1,pfpropiedades);
    fwrite(&Seg[i].Mass[0],sizeof(type_real),2,pfpropiedades); 
    fwrite(&Seg[i].Vnodos[0],sizeof(type_real),6,pfpropiedades); 
    fwrite(&Seg[i].razon,sizeof(type_real),1,pfpropiedades); 
    fwrite(&Seg[i].len,sizeof(type_real),1,pfpropiedades);
    fwrite(&Seg[i].elong,sizeof(type_real),1,pfpropiedades);
    fwrite(&Seg[i].rms,sizeof(type_real),1,pfpropiedades);

    fwrite(&Seg[i].size,sizeof(type_int),1,pfout);

    k = 0;
    for(j=Seg[i].start;j<Seg[i].start+Seg[i].size;j++)
    {
      fwrite(&Seg_centr[j].Pos[0],sizeof(type_real),1,pfout); 
      fwrite(&Seg_centr[j].Pos[1],sizeof(type_real),1,pfout); 
      fwrite(&Seg_centr[j].Pos[2],sizeof(type_real),1,pfout); 

      if(Seg[i].k_neigh[k].size()<K_NEIGH_SIZE)
      {
        assert(Seg[i].k_neigh[k].size()>0);
        fprintf(stdout,"%u %u tiene menos de %lu vecinos %u\n",i,j,Seg[i].k_neigh[k].size(),K_NEIGH_SIZE);
      }else if(Seg[i].k_neigh[k].size()>K_NEIGH_SIZE) {
        realloc_set(Seg[i].k_neigh[k]);
      }

      len = mass = 0.0f;
      for(std::set<std::pair<type_real,type_int> >::iterator it = Seg[i].k_neigh[k].begin(); \
        it != Seg[i].k_neigh[k].end(); ++it)
      {
        rs = (*it);
        len = rs.first; // se guarda la distancia al mas lejano
        #ifdef PARTICLES
        mass += cp.Mpart;
        #else
        mass += Gr[rs.second].NumPart*cp.Mpart;
        #endif  
      }
  
      len  /= 1000.;
      dens = (3.f*mass)/(4.*len*len*len); // Masa [10^10 Msol / h] / Volumen Mpc^3

      fwrite(&dens,sizeof(type_real),1,pfout); 

      Seg[i].k_neigh[k].empty();
      k++;
    }

    omp_destroy_lock(&(lock[i]));
  }

  fclose(pfout);
  fclose(pfpropiedades);
  free(lock);

  return;
}
