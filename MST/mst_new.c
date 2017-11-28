#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "voronoi.h"
#include "mst_kruskal.h"
#include "colores.h"

#ifdef CALCULA_MEDIA
  #include "grid.h"
#endif

int  NumPartCut;
std::vector<std::pair<int,int> > mass_orden;
#ifdef BRANCH_SURVIVE
  int N_part_survive;
#endif

void Write_Segments(int *Padre, int *Rank, type_real *fof);
#ifdef CALCULA_MEDIA
 void calc_media(type_real xc, type_real yc, type_real zc, \
  type_real vx, type_real vy, type_real vz, int *numpart, type_real *vel_media, \
  type_real *rcil2, type_real rlong_2, int ncil);
 bool point_inside(type_real dot, type_real rlong_2);
#endif

int main(int argc, char **argv)
{
  int    i;
  int    *Padre, *Rank;
  double start,end;
  double MassCut;
  std::vector<std::vector<int> > adjacency_list;
  std::vector<std::pair<type_real,std::pair<int,int> > > edges;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  /*Lee archivos de la simulacion*/
  read_gadget();

  MassCut = atof(argv[2]);
  NumPartCut = (int)(MassCut/cp.Mpart) ;  // Masa de la partícula [10^10 Msol / h]
  #ifdef BRANCH_SURVIVE
  N_part_survive = NumPartCut;
  #endif
  
  GREEN("********** Important *************\n");
  sprintf(message,"Mpart %g\n",cp.Mpart*1.e10);RED(message);
  sprintf(message,"Mass cut %g -> Nodos Npart %d Mass %g\n", \
  MassCut*1.e10,NumPartCut,NumPartCut*cp.Mpart*1e10);RED(message);
  GREEN("**********************************\n");
  fflush(stdout);

  select_particles_fof(fof[0]);

  read_grup_fof(fof[1]);

  Voronoi_Grupos(fof[0],edges);

  fprintf(stdout,"%d NumEdges\n",(int)edges.size());
  fflush(stdout);

  Padre = (int *) malloc(ngrup*sizeof(int));
  Rank =  (int *) malloc(ngrup*sizeof(int));

  for(i=0;i<ngrup;i++)
  {
    Padre[i] = i;
    Rank[i] = 0;
    adjacency_list.push_back(std::vector<int>());
  }

  Kruskal(Padre,Rank,edges,adjacency_list);

  Podado(4,adjacency_list);  

  #ifdef BRANCH_SURVIVE
  fprintf(stdout,"Sobreviven en ramas los nodos con %d\n",N_part_survive);
  fflush(stdout);
  #endif

  for(i=0;i<ngrup;i++)
  {
    Padre[i] = -1;
    Rank[i]  = 0;

    if(adjacency_list[i].size()>0)
      mass_orden.push_back(std::make_pair(Gr[i].NumPart,i));
  }

  sort(mass_orden.begin(),mass_orden.end());
  
  for(i=(int)mass_orden.size()-1;i>=0;i--)
  {
    int k = mass_orden[i].second;

    if(Padre[k]==-1)
      DLU(k,Padre[k],adjacency_list,Padre,Rank);
  }

  fprintf(stdout,"Escribe\n");
  fflush(stdout);

  Write_Segments(Padre,Rank,fof);

  free(Gr);
  free(Padre);
  free(Rank);
  adjacency_list.clear();
  #ifdef CALCULA_MEDIA
  free(P);
  #endif

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}

void Write_Segments(int *Padre, int *Rank, type_real *fof)
{
  char filename[200];
  FILE *pfout, *pfpropiedades;
  int i,j,k,id;
  type_real dx,dy,dz;
  type_real dux,duy,duz;
  type_real rms,elong;
  type_real r,lenr;
  std::vector<int> aux;
  std::vector<std::vector<int> > segmentos;
  #ifdef CALCULA_MEDIA
  int ncil;
  type_real xc,yc,zc;
  type_real rbin,racum;
  int bin,lbin;
  int *npart;
  type_real *Vmedia;
  type_real *rcil2;  // En Kpc
  type_real rlong = 500.; // En Kpc  
  FILE **pfvel;
  int   *c;

  ////////////////////////////////////////////////////////////////////////
  sprintf(filename,"Type in a number of bin\n"); GREEN(filename);
	scanf("%d",&ncil);
  sprintf(filename,"The number you typed was %d\n",ncil); GREEN(filename);
  fflush(stdout);
  ////////////////////////////////////////////////////////////////////////

  rcil2 = (type_real *) malloc(ncil*sizeof(type_real));

  for(i=0;i<ncil;i++)
  {
    r  = (type_real)(ncil-i)*1000.;
    rcil2[i] = r*r;
  }

  BLUE("********** Importante ***********\n");
  for(i=0;i<ncil;i++)
    {sprintf(filename,"Radio %d   %f Mpc\n",i,sqrt(rcil2[i])/1000.); RED(filename);}
  sprintf(filename,"Separacion %f Mpc\n",rlong/1000.);RED(filename);
  BLUE("**********************************\n");

  #endif

  j = 0;

  while(!mass_orden.empty())
  {

    i = mass_orden.back().second;

    //if(Padre[i]>=0)
    //if(Rank[i]==1 && Padre[i]>=0)
    if(Padre[i]>=0)
    {
      aux.push_back(i);
      id = Padre[i];
      Rank[i] *= -1;

      while(id>=0)
      {
        if(Rank[id]<0)
        {
          aux.push_back(id);
          break;
        }

        //if(Rank[id]>2)
        if(Gr[id].NumPart>NumPartCut && Padre[id]!=-1)
        {
          aux.push_back(id);

          segmentos.push_back(aux);
  
          aux.clear();
          j++;
        }

        aux.push_back(id);
        Rank[id] *= -1;
        id = Padre[id];               
      }

      segmentos.push_back(aux);

      aux.clear();
      j++;
    }

    mass_orden.pop_back();
  }

  assert((int)segmentos.size()==j);

  #ifdef MCRITIC
  sprintf(filename,"%.2d_%.4d_segmentos_cut_%.2f_%.2f_%.2f.bin",snap.num,NumPartCut,m_critica,fof[0],fof[1]);
  #else
  sprintf(filename,"%.2d_%.4d_segmentos_%.2f_%.2f.bin",snap.num,NumPartCut,fof[0],fof[1]);
  #endif
  pfout=fopen(filename,"w");
  fwrite(&j,sizeof(int),1,pfout);

  #ifdef MCRITIC
  sprintf(filename,"%.2d_%.4d_propiedades_cut_%.2f_%.2f_%.2f.bin",snap.num,NumPartCut,m_critica,fof[0],fof[1]);
  #else
  sprintf(filename,"%.2d_%.4d_propiedades_%.2f_%.2f.bin",snap.num,NumPartCut,fof[0],fof[1]);
  #endif

  pfpropiedades=fopen(filename,"w");
  fwrite(&j,sizeof(int),1,pfpropiedades);

  #ifdef SORT_DERECHA

  BLUE("********** Importante ***********\n");
  sprintf(filename,"Ordena el nodo de la derecha es más grande\n");RED(filename);
  BLUE("*********************************\n");

  #endif

  for(i=0;i<j;i++)
  {
    #ifdef SORT_DERECHA
    if(Gr[segmentos[i][0]].NumPart>Gr[segmentos[i].back()].NumPart)
      reverse(segmentos[i].begin(),segmentos[i].end());
    #endif

    aux = segmentos[i];

    #ifndef CALCULA_MEDIA
    segmentos[i].clear();
    #endif

    id = aux.back();

    k = 0;
    if(Gr[aux[0]].NumPart>NumPartCut) k++;
    if(Gr[id].NumPart>NumPartCut)     k++;

    fwrite(&k,sizeof(int),1,pfpropiedades);

    k = (int)aux.size();
    fwrite(&k,sizeof(int),1,pfout);
    fwrite(&k,sizeof(int),1,pfpropiedades);

    r = (float)Gr[aux[0]].NumPart/(float)Gr[id].NumPart;
    fwrite(&r,sizeof(float),1,pfpropiedades); // ESCRIBE LA RAZON DE MASAS

    dux = Gr[id].Pos[0] - Gr[aux[0]].Pos[0];
    duy = Gr[id].Pos[1] - Gr[aux[0]].Pos[1];
    duz = Gr[id].Pos[2] - Gr[aux[0]].Pos[2];

    #ifdef PERIODIC
    dux = dux >= cp.lbox*0.5 ? dux-cp.lbox : dux;
    dux = dux < -cp.lbox*0.5 ? dux+cp.lbox : dux;

    duy = duy >= cp.lbox*0.5 ? duy-cp.lbox : duy;
    duy = duy < -cp.lbox*0.5 ? duy+cp.lbox : duy;

    duz = duz >= cp.lbox*0.5 ? duz-cp.lbox : duz;
    duz = duz < -cp.lbox*0.5 ? duz+cp.lbox : duz;
    #endif

    elong = dux*dux+duy*duy+duz*duz;

    lenr = rms = 0.0;

    fwrite(&aux[0],sizeof(int),1,pfout);

    for(k=1;k<(int)aux.size();k++)
    {
      fwrite(&aux[k],sizeof(int),1,pfout);

      dx = Gr[aux[k]].Pos[0] - Gr[aux[k-1]].Pos[0];
      dy = Gr[aux[k]].Pos[1] - Gr[aux[k-1]].Pos[1];
      dz = Gr[aux[k]].Pos[2] - Gr[aux[k-1]].Pos[2];

      #ifdef PERIODIC
      dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
      dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

      dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
      dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

      dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
      dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
      #endif
      
      r = sqrt(dx*dx+dy*dy+dz*dz);

      lenr += r;

      if(k==(int)aux.size()-1) continue;

      dx = Gr[aux[k]].Pos[0] - Gr[aux[0]].Pos[0];
      dy = Gr[aux[k]].Pos[1] - Gr[aux[0]].Pos[1];
      dz = Gr[aux[k]].Pos[2] - Gr[aux[0]].Pos[2];

      #ifdef PERIODIC
      dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
      dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

      dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
      dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

      dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
      dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
      #endif
          
      r = pow(dy*duz-dz*duy,2);
      r += pow(dz*dux-dx*duz,2);
      r += pow(dx*duy-dy*dux,2);
      r /= elong;
      
      rms += r;

    }

    r = sqrt(elong)/lenr;

    assert(r<1.+1e-06);

    k = (int)aux.size();
    rms /= (float)k;
    rms = sqrt(rms);

    fwrite(&lenr,sizeof(float),1,pfpropiedades);
    fwrite(&r,sizeof(float),1,pfpropiedades);
    fwrite(&rms,sizeof(float),1,pfpropiedades);

    aux.clear();
  }

  fclose(pfout);
  fclose(pfpropiedades);

  #ifdef CALCULA_MEDIA

  //////////////////////////////////////////////////
  fprintf(stdout,"Build grid\n");

  grid.step = 0; /// IMPORTANTE
  grid.nobj = cp.npart;
  r = sqrt(rcil2[0]);
 
  if(rlong>r){
    grid.ngrid = (int)(cp.lbox/rlong);
  }else{ 
    grid.ngrid = (int)(cp.lbox/r);
  }

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }

  grid_init();
  grid_build();

  fprintf(stdout,"End build grid\n");

  //////////////////////////////////////////////////
  pfvel = (FILE **) malloc(NTHREADS*sizeof(FILE));
  c = (int *) malloc(NTHREADS*sizeof(int));

  for(i=0;i<NTHREADS;i++)
  {
    c[i] = 0;
    #ifdef MCRITIC
    sprintf(filename,"%.2d_%.4d_vmedia_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NumPartCut,m_critica,fof[0],fof[1],i);
    #else
    sprintf(filename,"%.2d_%.4d_vmedia_%.2f_%.2f.%.2d.bin",snap.num,NumPartCut,fof[0],fof[1],i);
    #endif
    pfvel[i]=fopen(filename,"w");
    fwrite(&c[i],sizeof(int),1,pfvel[i]);        
    fwrite(&ncil,sizeof(int),1,pfvel[i]);
  }
  ////////////////////////////////////////////////
  
  elong = rlong;

  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(static) default(none) private(i,k,aux,dx,dy,dz,r, \
  npart,Vmedia,dux,duy,duz,racum,bin,lbin,rbin,xc,yc,zc) \
  shared(j,pfvel,segmentos,Gr,cp,rlong,elong,rcil2,c,ncil,stdout)
  for(i=0;i<j;i++)
  {
    aux = segmentos[i];
    segmentos[i].clear();

    int Tid = omp_get_thread_num();
    npart = (int *) calloc(ncil,sizeof(int));
    Vmedia = (type_real *) calloc(3*ncil,sizeof(type_real));
    racum = 0.0;

    for(k=1;k<(int)aux.size();k++)
    {

      dx = Gr[aux[k]].Pos[0] - Gr[aux[k-1]].Pos[0];
      dy = Gr[aux[k]].Pos[1] - Gr[aux[k-1]].Pos[1];
      dz = Gr[aux[k]].Pos[2] - Gr[aux[k-1]].Pos[2];

      #ifdef PERIODIC
      dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
      dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

      dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
      dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

      dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
      dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
      #endif
      
      r = sqrt(dx*dx+dy*dy+dz*dz);

      dx /= r;
      dy /= r;
      dz /= r;

      if(k==1)
      {

        xc = Gr[aux[0]].Pos[0];
        yc = Gr[aux[0]].Pos[1];
        zc = Gr[aux[0]].Pos[2]; 

        calc_media(xc,yc,zc,dx,dy,dz,npart,Vmedia,rcil2,elong,ncil);            

      }else{

        xc = Gr[aux[k-1]].Pos[0];
        yc = Gr[aux[k-1]].Pos[1];
        zc = Gr[aux[k-1]].Pos[2];

        rbin = rlong-racum;
        xc += rbin*dx;
        yc += rbin*dy;
        zc += rbin*dz;   

        #ifdef PERIODIC
        xc = xc >= cp.lbox ? xc-cp.lbox : xc;
        xc = xc <      0.0 ? xc+cp.lbox : xc;
      
        yc = yc >= cp.lbox ? yc-cp.lbox : yc;
        yc = yc <      0.0 ? yc+cp.lbox : yc;
 
        zc = zc >= cp.lbox ? zc-cp.lbox : zc;
        zc = zc <      0.0 ? zc+cp.lbox : zc;
        #endif

        calc_media(xc,yc,zc,dx,dy,dz,npart,Vmedia,rcil2,elong,ncil);

        dux = Gr[aux[k]].Pos[0]-xc;
        duy = Gr[aux[k]].Pos[1]-yc;
        duz = Gr[aux[k]].Pos[2]-zc;

        #ifdef PERIODIC
        dux = dux >= 0.5*cp.lbox ? dux-cp.lbox : dux;
        dux = dux < -0.5*cp.lbox ? dux+cp.lbox : dux;
      
        duy = duy >= 0.5*cp.lbox ? duy-cp.lbox : duy;
        duy = duy < -0.5*cp.lbox ? duy+cp.lbox : duy;
 
        duz = duz >= 0.5*cp.lbox ? duz-cp.lbox : duz;
        duz = duz < -0.5*cp.lbox ? duz+cp.lbox : duz;
        #endif

        r = sqrt(dux*dux+duy*duy+duz*duz);

      }

      lbin = (int)(r/rlong);

      for(bin=0;bin<lbin;bin++)
      {
        xc += rlong*dx;
        yc += rlong*dy;
        zc += rlong*dz;   

        #ifdef PERIODIC
        xc = xc >= cp.lbox ? xc-cp.lbox : xc;
        xc = xc <      0.0 ? xc+cp.lbox : xc;
        
        yc = yc >= cp.lbox ? yc-cp.lbox : yc;
        yc = yc <      0.0 ? yc+cp.lbox : yc;
 
        zc = zc >= cp.lbox ? zc-cp.lbox : zc;
        zc = zc <      0.0 ? zc+cp.lbox : zc;
        #endif

        calc_media(xc,yc,zc,dx,dy,dz,npart,Vmedia,rcil2,elong,ncil);            
      }

      dux = xc-Gr[aux[k]].Pos[0];
      duy = yc-Gr[aux[k]].Pos[1];
      duz = zc-Gr[aux[k]].Pos[2];

      #ifdef PERIODIC
      dux = dux >= 0.5*cp.lbox ? dux-cp.lbox : dux;
      dux = dux < -0.5*cp.lbox ? dux+cp.lbox : dux;
      
      duy = duy >= 0.5*cp.lbox ? duy-cp.lbox : duy;
      duy = duy < -0.5*cp.lbox ? duy+cp.lbox : duy;
 
      duz = duz >= 0.5*cp.lbox ? duz-cp.lbox : duz;
      duz = duz < -0.5*cp.lbox ? duz+cp.lbox : duz;
      #endif
      
      racum = sqrt(dux*dux+duy*duy+duz*duz);
      
    }

    calc_media(Gr[aux.back()].Pos[0],Gr[aux.back()].Pos[1],Gr[aux.back()].Pos[2], \
    dx,dy,dz,npart,Vmedia,rcil2,elong,ncil);            

    ///////////////////////////////////////////////////////////////////////////////////

    dx = Gr[aux[1]].Pos[0] - Gr[aux[0]].Pos[0];
    dy = Gr[aux[1]].Pos[1] - Gr[aux[0]].Pos[1];
    dz = Gr[aux[1]].Pos[2] - Gr[aux[0]].Pos[2];

    #ifdef PERIODIC
    if(dx> 0.5*cp.lbox) dx -= cp.lbox;
    if(dx<-0.5*cp.lbox) dx += cp.lbox;
    if(dy> 0.5*cp.lbox) dy -= cp.lbox;
    if(dy<-0.5*cp.lbox) dy += cp.lbox;
    if(dz> 0.5*cp.lbox) dz -= cp.lbox;
    if(dz<-0.5*cp.lbox) dz += cp.lbox;
    #endif

    r = sqrt(dx*dx+dy*dy+dz*dz);

    dx /= r;
    dy /= r;
    dz /= r;

    xc = Gr[aux[0]].Pos[0]+racum*dx;
    yc = Gr[aux[0]].Pos[1]+racum*dy;
    zc = Gr[aux[0]].Pos[2]+racum*dz;

    #ifdef PERIODIC
    xc = xc >= cp.lbox ? xc-cp.lbox : xc;
    xc = xc <      0.0 ? xc+cp.lbox : xc;
    
    yc = yc >= cp.lbox ? yc-cp.lbox : yc;
    yc = yc <      0.0 ? yc+cp.lbox : yc;
 
    zc = zc >= cp.lbox ? zc-cp.lbox : zc;
    zc = zc <      0.0 ? zc+cp.lbox : zc;
    #endif

    calc_media(xc,yc,zc,dx,dy,dz,npart,Vmedia,rcil2,elong,ncil);             
    ///////////////////////////////////////////////////////////////////////////////////
    
    for(bin=0;bin<ncil;bin++)
    {
      for(k=0;k<3;k++)
      {
        if(npart[bin]==0)          
          fprintf(stdout,"Tid %d Fil %d BREAK!!! npart[bin]==0!!!\n",Tid,i);  
        Vmedia[bin+ncil*k] /= (type_real)npart[bin];
      }
    }

    fwrite(&i,sizeof(int),1,pfvel[Tid]);
    fwrite(&Vmedia[0],sizeof(type_real),3*ncil,pfvel[Tid]);
    c[Tid]++;

    free(npart);
    free(Vmedia);
    aux.clear();
  } // cierra el for

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfvel[i]);
    fwrite(&c[i],sizeof(int),1,pfvel[i]);
    fclose(pfvel[i]);
  }

  free(c);
  grid_free();
  free(rcil2);

  #endif // CALCULA MEDIA

  while(!segmentos.empty())
    segmentos.pop_back();

  fprintf(stdout,"Segmentos %d\n",j);

  return;

}

#ifdef CALCULA_MEDIA

void calc_media(type_real xc, type_real yc, type_real zc, \
type_real vx, type_real vy, type_real vz, int *numpart, type_real *vel_media, \
type_real *rcil2, type_real rlong_2, int ncil)
{
  int i,j,k;
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  type_real Posprima[3];
  type_real lbox,fac,lbox2;
  type_real dot,dis;
  int ngrid;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  ixc  = (int)(xc*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)(yc*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)(zc*fac);
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
          Posprima[0] = P[i].Pos[0] - xc;
          Posprima[1] = P[i].Pos[1] - yc;
          Posprima[2] = P[i].Pos[2] - zc;

          #ifdef PERIODIC
          if(Posprima[0] >  lbox2) Posprima[0] = Posprima[0] - lbox;
          if(Posprima[1] >  lbox2) Posprima[1] = Posprima[1] - lbox;
          if(Posprima[2] >  lbox2) Posprima[2] = Posprima[2] - lbox;
          if(Posprima[0] < -lbox2) Posprima[0] = Posprima[0] + lbox;
          if(Posprima[1] < -lbox2) Posprima[1] = Posprima[1] + lbox;
          if(Posprima[2] < -lbox2) Posprima[2] = Posprima[2] + lbox;
          #endif

          dot = Posprima[0]*vx+Posprima[1]*vy+Posprima[2]*vz;

          if(point_inside(dot,rlong_2))
          {        

            dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2]-dot*dot;        

            for(j=0;j<ncil;j++)
            {
              if(dis<rcil2[j])
              {

                for(k=0;k<3;k++)
                  vel_media[j+ncil*k] += P[i].Vel[k];

                numpart[j]++;

              }else{

                break;

              }
            }
          } // cierra dis


          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

bool point_inside(type_real dot, type_real rlong_2)
{

  if(fabs(dot)>rlong_2){

      return false;

  }else{

    return true;

  }

}
#endif
