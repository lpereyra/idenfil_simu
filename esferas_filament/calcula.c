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

struct gridst grid;

#ifdef SPH

  static type_real std_normal(type_real x_2)
  {
      return exp(-0.5*x_2);     //exp(-0.5*x*x);
  }

#else

  static void realloc_set(std::set<std::pair<type_real,type_int> > &myset)
  {
    std::set<std::pair<type_real,type_int> >::iterator it = myset.begin(); 
    std::advance (it,K_NEIGH_SIZE);
    myset.erase(it,myset.end());
  
    return;
  }

#endif


#ifdef SPH
  static void calc_search_fil(const type_real * __restrict__ Pos_cent, type_real &dens, type_real &vol, \
                              const type_int rint, const type_real r_2) 
#else
  static void calc_search_fil(const type_real * __restrict__ Pos_cent, std::set<std::pair<type_real,type_int> > &myset, \
                              const type_int rint, const type_real r_2) 
#endif
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
  ixci = ixc - rint;
  ixcf = ixc + rint;
  iyc  = (long)(Pos_cent[1]*fac);
  iyci = iyc - rint;
  iycf = iyc + rint;
  izc  = (long)(Pos_cent[2]*fac);
  izci = izc - rint;
  izcf = izc + rint;

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
          #ifdef PARTICLES
            Posprima[0] = P[i].Pos[0] - Pos_cent[0];
            Posprima[1] = P[i].Pos[1] - Pos_cent[1];
            Posprima[2] = P[i].Pos[2] - Pos_cent[2];
          #else
            Posprima[0] = Gr[i].Pos[0] - Pos_cent[0];
            Posprima[1] = Gr[i].Pos[1] - Pos_cent[1];
            Posprima[2] = Gr[i].Pos[2] - Pos_cent[2];
          #endif

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
            #ifdef SPH

              #ifdef PARTICLES
                vol  += cp.Mpart;
                dens += cp.Mpart*std_normal(dis/H_SOFT_2);
              #else
                vol  += Gr[i].NumPart*cp.Mpart;
                dens += Gr[i].NumPart*cp.Mpart*std_normal(dis/H_SOFT_2);
              #endif  

            #else

              if(myset.size()>SET_MAX_SIZE)
                realloc_set(myset);

              myset.insert(std::make_pair(dis,i)); //dis_2

            #endif
          }// cierra el if

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

extern void propiedades(type_int NNN, type_real *fof, const type_real r_aux, const type_int rint, type_int &flag)
{
  type_int i,j;
  const type_real rbus_2 = ((type_real)rint*r_aux)*((type_real)rint*r_aux);
  #ifdef SPH
  const type_real norm_cte = pow(H_SOFT/1000.,3)*0.5*M_SQRT1_2*M_2_SQRTPI; // (H_SOFT/1000.^3.0) * (0.5*1/sqrt(2)*2/sqrt(pi) || sqrt(2*M_PI))
  #else
  const type_real norm_cte = (4.*M_PI/3.); // (4*pi/3)
  #endif

  fprintf(stdout,"iteracion %d faltan %d rbus %f\n",rint,flag,sqrt(rbus_2));

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) private(i,j) \
  shared(cp,Seg,Seg_centr,Gr,stdout) reduction(-:flag)
  //for(i=0;i<cp.nseg;i++)
  for(i=0;i<cp.ncent_seg;i++)
  {
    //if(i%100000==0) 
    if(i%1000000==0) 
    {
      //fprintf(stdout,"%u %.4f\n",i,(type_real)i/(type_real)cp.nseg);
      fprintf(stdout,"%u %.4f\n",i,(type_real)i/(type_real)cp.ncent_seg);
      fflush(stdout);
    }

    //for(j=0;j<Seg[i].size;j++)
    //{
 
      j = Seg_centr[i].id;
      type_real vol = 0.0f;

      #ifdef SPH

        flag -= 1;
        //Seg[i].dens[j] =  0.0;    
        //calc_search_fil(Seg_centr[Seg[i].start+j].Pos,Seg[i].dens[j],vol,rint,rbus_2);
        Seg[Seg_centr[i].save].dens[j] =  0.0;    
        calc_search_fil(Seg_centr[i].Pos,Seg[Seg_centr[i].save].dens[j],vol,rint,rbus_2);

        vol = vol*norm_cte;

        //Seg[i].dens[j] /= vol; // Masa [10^10 Msol / h] / Volumen Mpc^3
        Seg[Seg_centr[i].save].dens[j] /= vol; // Masa [10^10 Msol / h] / Volumen Mpc^3

      #else

        if(Seg[Seg_centr[i].save].dens[j]>0.0f) continue;

        std::set<std::pair<type_real,type_int> > k_neigh;

        calc_search_fil(Seg_centr[i].Pos,k_neigh,rint,rbus_2);

        if(k_neigh.size()>=K_NEIGH_SIZE) 
        {
          flag -= 1;

          //Seg[i].dens[j] = 0.0f;          
          Seg[Seg_centr[i].save].dens[j] = 0.0f;

          realloc_set(k_neigh);

          for(std::set<std::pair<type_real,type_int> >::iterator it = k_neigh.begin(); \
            it != k_neigh.end(); ++it)
          {
            std::pair<type_real, type_int> rs = (*it);

            vol = rs.first; // se guarda la distancia al mas lejano
            #ifdef PARTICLES
              //Seg[i].dens[j] += cp.Mpart;
              Seg[Seg_centr[i].save].dens[j] += cp.Mpart;
            #else
              //Seg[i].dens[j] += cp.Mpart*(type_real)Gr[rs.second].NumPart;
              Seg[Seg_centr[i].save].dens[j] += cp.Mpart*(type_real)Gr[rs.second].NumPart;
            #endif  
          }

          vol /= pow(1000.,2); // (1000.^2.0)
          vol = sqrt(vol)*vol*norm_cte;

          //Seg[i].dens[j] /= vol; // Masa [10^10 Msol / h] / Volumen Mpc^3
          Seg[Seg_centr[i].save].dens[j] /= vol; // Masa [10^10 Msol / h] / Volumen Mpc^3
        }

        k_neigh.empty();

      #endif

    //}
  }

  BLUE("******************************\n");
  GREEN("******************************\n");
  fflush(stdout);
 
  return;
}

extern void set_name(char * name, const type_int NNN, const char * prefix)
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

  #ifdef PARTICLES
    sprintf(name,"%s_PARTICULAS",name);
  #else
    sprintf(name,"%s_HALOS",name);
  #endif

  #ifdef SPH
    sprintf(name,"%s_sph_%.2f",name,H_SOFT/1000.);
  #else
    sprintf(name,"%s_k_neigh_%.2d",name,K_NEIGH_SIZE);
  #endif

  sprintf(name,"%s_%.2f_%.2f.bin",name,fof[0],fof[1]);

  return;
}

