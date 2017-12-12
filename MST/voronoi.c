#include <assert.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <functional>
#include <algorithm>
#include <vector>

#include "cosmoparam.h"
#include "variables.h"
#include "grid.h"
#include "voronoi.h"
#include "voro++.hh"
#ifdef WRITE_WEIGHT
  #include "leesnap.h"
#endif

bool dist_segment(double fof_2, int idg, type_real x, type_real y, type_real z)
{
  int i;
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  type_real lbox,fac,lbox2;
  double dx,dy,dz,r;
  int ngrid;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  ixc  = (int)(x*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)(y*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)(z*fac);
  izci = izc - 1;
  izcf = izc + 1;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= ngrid ) ixcf = ngrid - 1;
  if( iycf >= ngrid ) iycf = ngrid - 1;
  if( izcf >= ngrid ) izcf = ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= ngrid) ix = ix - ngrid;
    if(ix < 0) ix = ix + ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= ngrid) iy = iy - ngrid;
      if(iy < 0) iy = iy + ngrid;
      #endif
  
      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= ngrid) iz = iz - ngrid;
        if(iz < 0) iz = iz + ngrid;
        #endif

        ibox = (ix * ngrid + iy) * ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {

          if(P[i].sub==idg)
          {
            dx = P[i].Pos[0] - x;
            dy = P[i].Pos[1] - y;
            dz = P[i].Pos[2] - z;

            #ifdef PERIODIC
            dx = dx >= lbox2 ? dx-lbox : dx;
            dx = dx < -lbox2 ? dx+lbox : dx;

            dy = dy >= lbox2 ? dy-lbox : dy;
            dy = dy < -lbox2 ? dy+lbox : dy;
  
            dz = dz >= lbox2 ? dz-lbox : dz;
            dz = dz < -lbox2 ? dz+lbox : dz;
            #endif

            r = dx*dx+dy*dy+dz*dz;

            if(r<fof_2)
              return true;
          }

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return false;
}

void Voronoi_Grupos(type_real fof, std::vector<std::pair<type_real,std::pair<int,int> > > &edges)
{

  int i, j, idv, id, Ngrid, Tid, itera;
  bool xbool, ybool, zbool;
  double x_min, y_min, z_min;	  
  double x_max, y_max, z_max;	  
  double dx,dy,dz,r,r0,r0_2,frac;
  std::vector<int>  vec;
  voro::voronoicell_neighbor cell;
  #ifdef WRITE_WEIGHT
    char filename[200];
    int  *c;
    FILE **pfpesos;
    pfpesos = (FILE **) malloc(NTHREADS*sizeof(FILE));
    c       = (int  *)  malloc(NTHREADS*sizeof(int));
  #endif

  Ngrid = (int)pow((float)ngrup/5.0,1.0/3.0);
  #ifdef PERIODIC
    xbool = ybool = zbool = true;
  #else
    xbool = ybool = zbool = false;
  #endif
  r0 = 1.0e-10; // uso un pequeño gap
  x_min = y_min = z_min = 0.0-r0;	  
  x_max = y_max = z_max = cp.lbox+r0;	  

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  r0 = fof*cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.;

  /// INICIA LOS PARAMETROS DE LA GRILLA ///
  grid.nobj = cp.npart;
  grid.ngrid = (int)(1.1*cp.lbox/r0);
  grid.step = 1;

  fprintf(stdout,"r0 %f lbox %f\n",r0,cp.lbox);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }
  //////////////////////////////////////////
  grid_init();
  grid_build();

  r0_2 = r0*r0; // Para usar el cuadrado
  r0 /= 100.;

  // Creo el container
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,Ngrid,Ngrid,Ngrid,xbool,ybool,zbool,8);

  for(i=0;i<ngrup; i++)
    con.put(i,Gr[i].Pos[0],Gr[i].Pos[1],Gr[i].Pos[2]);	  

  assert(ngrup==con.total_particles());

  voro::c_loop_all clo(con);
  std::vector<std::pair<int,int> > pares;
  std::vector<std::vector<std::pair<type_real,std::pair<int,int> > > > lados(NTHREADS);

  if(clo.start()) do if(con.compute_cell(cell,clo))
  {
    id = clo.pid();
    cell.neighbors(vec);

    for(j=0; j<(int)vec.size(); j++)
    {
      idv = vec[j];

      if(Gr[id].save!=Gr[idv].save) continue;

      if(Gr[id].id>Gr[idv].id)
        pares.push_back(std::make_pair(idv,id));

    }

    vec.clear();

  }while(clo.inc());

  con.clear();

  #ifdef WRITE_WEIGHT
    for(i=0;i<NTHREADS;i++)
    {
      c[i] = 0;
      sprintf(filename,"%.2d_pesos.%.2d.bin",snap.num,i);
      pfpesos[i] = fopen(filename,"w");
      fwrite(&c[i],sizeof(int),1,pfpesos[i]);    
    }


    #pragma omp parallel for num_threads(NTHREADS) schedule(dynamic) default(none) \
    private(i,Tid,id,idv,dx,dy,dz,r,itera,frac) \
    shared(r0,r0_2,vec,cp,Gr,grid,pares,lados,pfpesos,c,stdout)

  #else

    #pragma omp parallel for num_threads(NTHREADS) schedule(dynamic) default(none) \
    private(i,Tid,id,idv,dx,dy,dz,r,itera,frac) \
    shared(r0,r0_2,vec,cp,Gr,grid,pares,lados,stdout)

  #endif
  for(i=0;i<(int)pares.size();i++)
  {

    Tid = omp_get_thread_num();
    id  = pares[i].first;
    idv = pares[i].second;

    dx = Gr[idv].Pos[0] - Gr[id].Pos[0];
    dy = Gr[idv].Pos[1] - Gr[id].Pos[1];
    dz = Gr[idv].Pos[2] - Gr[id].Pos[2];

    #ifdef PERIODIC
    dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
    dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

    dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
    dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

    dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
    dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
    #endif

    r = sqrt(dx*dx+dy*dy+dz*dz);
    itera = (int)(r/r0);

    dx /= r;
    dy /= r;
    dz /= r;

    while(itera>0)
    {
      frac = (type_real)(itera*r0);
      if(!dist_segment(r0_2,Gr[id].save,Gr[id].Pos[0]+frac*dx,Gr[id].Pos[1]+frac*dy,Gr[id].Pos[2]+frac*dz))
        break;
      else  
        itera--; 
    }

    if(itera!=0) continue;

    #ifdef WRITE_WEIGHT
    fwrite(&id,sizeof(int),1,pfpesos[Tid]);
    fwrite(&idv,sizeof(int),1,pfpesos[Tid]);
    fwrite(&Gr[id].NumPart,sizeof(int),1,pfpesos[Tid]);
    fwrite(&Gr[idv].NumPart,sizeof(int),1,pfpesos[Tid]);
    fwrite(&r,sizeof(double),1,pfpesos[Tid]); // distancia
    #endif

    r = -((float)Gr[id].NumPart*(float)Gr[idv].NumPart)/(r*r);

    lados[Tid].push_back(std::make_pair((type_real)r,pares[i]));

    #ifdef WRITE_WEIGHT
    fwrite(&r,sizeof(double),1,pfpesos[Tid]); // peso
    c[Tid]++;
    #endif

  }

  pares.clear();

  #ifdef WRITE_WEIGHT
  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfpesos[i]);
    fwrite(&c[i],sizeof(int),1,pfpesos[i]);
    fclose(pfpesos[i]);
  }

  free(c);
  #endif

  for(i=0;i<NTHREADS;i++)
  {
    edges.insert(edges.end(),lados[i].begin(),lados[i].end());
    lados[i].clear();
  }

  grid_free();
  #ifndef CALCULA_MEDIA
  free(P);
  #endif
  
  return;
}
