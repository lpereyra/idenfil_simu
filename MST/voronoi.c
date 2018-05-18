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

  int i, idv, id, Tid;
  type_real r,r0,r0_2;
  type_real vdir[3];
  std::vector<int>  vec;
  std::vector<std::pair<int,int> > pares;
  std::vector<std::vector<std::pair<type_real,std::pair<int,int> > > > lados(NTHREADS);
  #ifdef WRITE_WEIGHT
    char filename[200];
    int  *c;
    FILE **pfpesos;
  #endif

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  ///////////////////////////////////////////

  init_edges(pares);

  ///////////////////////////////////////////
 
  //// INICIA LOS PARAMETROS DE LA GRILLA ////
  r0 = fof*cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.f;
  grid.step  = 1;
  grid.nobj  = cp.npart;
  grid.ngrid = (int)(1.1f*cp.lbox/r0);
  if(grid.ngrid > NGRIDMAX)
  {
    grid.ngrid = NGRIDMAX;
    fprintf(stdout,"Using NGRIDMAX = %lu, r0 = %f, lbox = %f\n",grid.ngrid,r0,cp.lbox);
  }else{
    fprintf(stdout,"Ngrid = %lu, r0 = %f, lbox = %f\n",grid.ngrid,r0,cp.lbox);
  }

  grid_init();
  grid_build();
  r0_2 = r0*r0; // Para usar el cuadrado
  ///////////////////////////////////////////

  #ifdef WRITE_WEIGHT
    c       = (int  *)  malloc(NTHREADS*sizeof(int));
    pfpesos = (FILE **) malloc(NTHREADS*sizeof(FILE));

    for(i=0;i<NTHREADS;i++)
    {
      c[i] = 0;
      sprintf(filename,"%.2d_pesos.%.2d.bin",snap.num,i);
      pfpesos[i] = fopen(filename,"w");
      fwrite(&c[i],sizeof(int),1,pfpesos[i]);    
    }
  #endif

  #ifdef WRITE_WEIGHT

    #pragma omp parallel for num_threads(NTHREADS) schedule(dynamic) default(none) \
    private(i,Tid,id,idv,r,vdir) \
    shared(fof,r0,r0_2,vec,cp,Gr,grid,pares,lados,pfpesos,c,stdout)

  #else

    #pragma omp parallel for num_threads(NTHREADS) schedule(dynamic) default(none) \
    private(i,Tid,id,idv,r,vdir) \
    shared(fof,r0,r0_2,vec,cp,Gr,grid,pares,lados,stdout)

  #endif
  for(i=0;i<(int)pares.size();i++)
  {

    Tid = omp_get_thread_num();
    id  = pares[i].first;
    idv = pares[i].second;

    vdir[0] = Gr[idv].Pos[0] - Gr[id].Pos[0];
    vdir[1] = Gr[idv].Pos[1] - Gr[id].Pos[1];
    vdir[2] = Gr[idv].Pos[2] - Gr[id].Pos[2];

    #ifdef PERIODIC
    vdir[0] = vdir[0] >= cp.lbox*0.5 ? vdir[0]-cp.lbox : vdir[0];
    vdir[0] = vdir[0] < -cp.lbox*0.5 ? vdir[0]+cp.lbox : vdir[0];

    vdir[1] = vdir[1] >= cp.lbox*0.5 ? vdir[1]-cp.lbox : vdir[1];
    vdir[1] = vdir[1] < -cp.lbox*0.5 ? vdir[1]+cp.lbox : vdir[1];

    vdir[2] = vdir[2] >= cp.lbox*0.5 ? vdir[2]-cp.lbox : vdir[2];
    vdir[2] = vdir[2] < -cp.lbox*0.5 ? vdir[2]+cp.lbox : vdir[2];
    #endif

    r = sqrt(vdir[0]*vdir[0]+vdir[1]*vdir[1]+vdir[2]*vdir[2]);

    vdir[0] *= (1./r);
    vdir[1] *= (1./r);
    vdir[2] *= (1./r);

    if(inside_fof(fof,id,r0,r0_2,r,vdir)==0) continue;

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

void init_edges(std::vector<std::pair<int,int> > &pares)
{

  int i, Ngrid;
  const double gap = 1.0e-10; // uso un pequeño gap
  double x_min, y_min, z_min;	  
  double x_max, y_max, z_max;	  
  bool xbool, ybool, zbool;
  voro::voronoicell_neighbor cell;

  Ngrid = (type_int)pow((float)cp.ngrup/5.0,1.0/3.0);
  #ifdef PERIODIC
    xbool = ybool = zbool = true;
  #else
    xbool = ybool = zbool = false;
  #endif

  x_min = y_min = z_min = 0.0-gap;	  
  x_max = y_max = z_max = cp.lbox+gap;	  

  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,Ngrid,Ngrid,Ngrid,xbool,ybool,zbool,8);

  for(i=0;i<cp.ngrup; i++)
  {
    con.put(i,Gr[i].Pos[0],Gr[i].Pos[1],Gr[i].Pos[2]);	  
  }
  assert(cp.ngrup==con.total_particles());

  voro::c_loop_all clo(con);

  if(clo.start()) do if(con.compute_cell(cell,clo))
  {
    std::vector<int>  vec;
    int id = clo.pid();
    cell.neighbors(vec);

    for(i=0; i<(int)vec.size(); i++)
    {
      int idv = vec[i];
    
      #ifndef PERIODIC
      if(idv<0) continue;
      #endif

      if(Gr[id].save != Gr[idv].save) continue;

      if(Gr[id].id > Gr[idv].id)
        pares.push_back(std::make_pair(idv,id));

    }

    vec.clear();

  }while(clo.inc());

  con.clear();

  return;
}

int inside_fof(const type_real fof, const type_int id, const type_real rcil, const type_real rcil2, const type_real r, const type_real vdir[])
{   

   type_real rsep   = rcil; // rcil == rsep -> rsph = 2**0.5 * rcil
   type_int  itera  = (int)(r/rsep);
   type_real rsph, Pos_cent[3];

   Pos_cent[0] = Gr[id].Pos[0];
   Pos_cent[1] = Gr[id].Pos[1];
   Pos_cent[2] = Gr[id].Pos[2];

   if(itera==0)
   {     
      Pos_cent[0] += (0.5*r)*vdir[0];
      Pos_cent[1] += (0.5*r)*vdir[1];
      Pos_cent[2] += (0.5*r)*vdir[2];

      #ifdef PERIODIC
      Pos_cent[0] = Pos_cent[0]<0 ? cp.lbox+(type_real)fmod(Pos_cent[0],cp.lbox) : (type_real)fmod(Pos_cent[0],cp.lbox);
      Pos_cent[1] = Pos_cent[1]<0 ? cp.lbox+(type_real)fmod(Pos_cent[1],cp.lbox) : (type_real)fmod(Pos_cent[1],cp.lbox);
      Pos_cent[2] = Pos_cent[2]<0 ? cp.lbox+(type_real)fmod(Pos_cent[2],cp.lbox) : (type_real)fmod(Pos_cent[2],cp.lbox);
      #endif

      rsph = sqrt(r*r+rcil2);

      if(!control(fof,0.5*r,rcil2,rsph,Gr[id].save,vdir,Pos_cent))
        return 0; // primer punto evaluado
   }

   rsph = sqrt(2.f)*rcil;
   Pos_cent[0] += (0.5f*rsep)*vdir[0];
   Pos_cent[1] += (0.5f*rsep)*vdir[1];
   Pos_cent[2] += (0.5f*rsep)*vdir[2];

   #ifdef PERIODIC
   Pos_cent[0] = Pos_cent[0]<0 ? cp.lbox+(type_real)fmod(Pos_cent[0],cp.lbox) : (type_real)fmod(Pos_cent[0],cp.lbox);
   Pos_cent[1] = Pos_cent[1]<0 ? cp.lbox+(type_real)fmod(Pos_cent[1],cp.lbox) : (type_real)fmod(Pos_cent[1],cp.lbox);
   Pos_cent[2] = Pos_cent[2]<0 ? cp.lbox+(type_real)fmod(Pos_cent[2],cp.lbox) : (type_real)fmod(Pos_cent[2],cp.lbox);
   #endif

   if(!control(fof,0.5f*rsep,rcil2,rsph,Gr[id].save,vdir,Pos_cent))
     return 0; // primer punto evaluado

   for(type_int c=1;c<itera;c++)
   {

     Pos_cent[0] += rsep*vdir[0];
     Pos_cent[1] += rsep*vdir[1];
     Pos_cent[2] += rsep*vdir[2];
     #ifdef PERIODIC
     Pos_cent[0] = Pos_cent[0]<0 ? cp.lbox+(type_real)fmod(Pos_cent[0],cp.lbox) : (type_real)fmod(Pos_cent[0],cp.lbox);
     Pos_cent[1] = Pos_cent[1]<0 ? cp.lbox+(type_real)fmod(Pos_cent[1],cp.lbox) : (type_real)fmod(Pos_cent[1],cp.lbox);
     Pos_cent[2] = Pos_cent[2]<0 ? cp.lbox+(type_real)fmod(Pos_cent[2],cp.lbox) : (type_real)fmod(Pos_cent[2],cp.lbox);
     #endif

     if(!control(fof,0.5*rsep,rcil2,rsph,Gr[id].save,vdir,Pos_cent))
       return 0; // primer punto evaluado

   }

   Pos_cent[0] += (0.5*rsep)*vdir[0];
   Pos_cent[1] += (0.5*rsep)*vdir[1];
   Pos_cent[2] += (0.5*rsep)*vdir[2];

   rsep = r - (type_real)itera*rsep; // cilindro pequeño centrado en la diferencia
   rsph = sqrt(rsep*rsep+rcil2);

   Pos_cent[0] += (0.5*rsep)*vdir[0];
   Pos_cent[1] += (0.5*rsep)*vdir[1];
   Pos_cent[2] += (0.5*rsep)*vdir[2];

   #ifdef PERIODIC
   Pos_cent[0] = Pos_cent[0]<0 ? cp.lbox+(type_real)fmod(Pos_cent[0],cp.lbox) : (type_real)fmod(Pos_cent[0],cp.lbox);
   Pos_cent[1] = Pos_cent[1]<0 ? cp.lbox+(type_real)fmod(Pos_cent[1],cp.lbox) : (type_real)fmod(Pos_cent[1],cp.lbox);
   Pos_cent[2] = Pos_cent[2]<0 ? cp.lbox+(type_real)fmod(Pos_cent[2],cp.lbox) : (type_real)fmod(Pos_cent[2],cp.lbox);
   #endif
  
   if(!control(fof,0.5*rsep,rcil2,rsph,Gr[id].save,vdir,Pos_cent))
     return 0; // ultimo punto evaluado

   return 1; 
}

inline long igrid(const long ix, const long iy, const long iz, const long ngrid) 
{
  return (ix*ngrid+ iy)*ngrid + iz;
}

int control(const type_real fof, const type_real rlong_2, const type_real rcil2, const type_real rsph, const int idg, const type_real versor[], const type_real Pos_cent[])
{
	
  int i, ilim;
  long ixc, iyc, izc, ibox;
  type_int  npart = 0;
  type_real Posprima[3];
  type_real dd, dis;

  ilim = (long)(rsph*(type_real)grid.ngrid*(1.f/cp.lbox))+1;
  ixc  = (long)(Pos_cent[0]*(type_real)grid.ngrid*(1.f/cp.lbox));
  iyc  = (long)(Pos_cent[1]*(type_real)grid.ngrid*(1.f/cp.lbox));
  izc  = (long)(Pos_cent[2]*(type_real)grid.ngrid*(1.f/cp.lbox));

  #ifndef PERIODIC
    for(long ixx = (ixc-ilim<0) ? 0 : ixc-ilim; ixx <= (ixc+ilim >= grid.ngrid) ? grid.ngrid-1 : ixc+ilim; ixx++)
  #else
    for(long ixx = ixc-ilim; ixx <= ixc+ilim; ixx++)
  #endif
  {
    #ifndef PERIODIC
      for(long iyy = (iyc-ilim<0) ? 0 : iyc-ilim; iyy <= (iyc+ilim >= grid.ngrid) ? grid.ngrid-1 : iyc+ilim; iyy++)
    #else
      for(long iyy = iyc-ilim ; iyy <= iyc+ilim ; iyy++)
    #endif
    {
      #ifndef PERIODIC
        for(long izz = (izc-ilim<0) ? 0 : izc-ilim; izz <= (izc+ilim >= grid.ngrid) ? grid.ngrid-1 : izc+ilim; izz++)
      #else
        for(long izz = izc-ilim ; izz <= izc+ilim ; izz++)
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

        i = grid.llirst[ibox];

        while(i != -1)
        {

          if(P[i].sub!=idg)
          {
            i = grid.ll[i];
            continue;  
          }
 
          Posprima[0] = P[i].Pos[0] - Pos_cent[0];
          Posprima[1] = P[i].Pos[1] - Pos_cent[1];
          Posprima[2] = P[i].Pos[2] - Pos_cent[2];

          #ifdef PERIODIC
          if(Posprima[0] >  cp.lbox*0.5f) Posprima[0] = Posprima[0] - cp.lbox;
          if(Posprima[1] >  cp.lbox*0.5f) Posprima[1] = Posprima[1] - cp.lbox;
          if(Posprima[2] >  cp.lbox*0.5f) Posprima[2] = Posprima[2] - cp.lbox;
          if(Posprima[0] < -cp.lbox*0.5f) Posprima[0] = Posprima[0] + cp.lbox;
          if(Posprima[1] < -cp.lbox*0.5f) Posprima[1] = Posprima[1] + cp.lbox;
          if(Posprima[2] < -cp.lbox*0.5f) Posprima[2] = Posprima[2] + cp.lbox;
          #endif

          dd = Posprima[0]*versor[0]+Posprima[1]*versor[1]+Posprima[2]*versor[2];

          if(fabs(dd)<rlong_2)
          {

            dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2] - dd*dd;

            if(dis<rcil2)
            { 
              npart++;
            }// cierra el if

          }

          i = grid.ll[i];

        } /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/  

  dd =  2.f*rlong_2*M_PI*rcil2/(type_real)npart/100.0f;
  dd = ((cp.lbox*cp.lbox*cp.lbox)/dd)/(type_real)cp.npart;

  if(dd>fof)
   return 1;
  else
   return 0;

}
