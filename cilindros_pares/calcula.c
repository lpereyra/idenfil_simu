#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "list.h"
#include "calcula.h"
#include "colores.h"
#include "leesnap.h"
#ifdef SAVEPART
  #include "bitmask.h"
#endif

void propiedades(int NNN, type_real *fof)
{
  char filename[200], name[200];
  int i, j, l, idim, Tid, *c;
  int totfil = 0;
  int binsper;
  type_real *rcil2;                 // En Kpc                    
  type_real rsep;                   // En Kpc
  type_real rlong;                  // En Kpc
  type_real rlong_2;                // En Kpc
  type_real RMAX, RMIN, LMAX;
  type_real r;
  FILE **pfdens;
  #ifdef CALCULA_MEDIA
  FILE **pfvel;
  #endif
  #ifdef SAVEPART
  FILE **pfpart;
  #endif
  #ifdef EXTEND
  FILE **pfextend;
  #endif
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
  binsper = ncil+1;

  //////////////////////////////////////////////////
  fprintf(stdout,"Build grid\n");

  grid.nobj = cp.npart;
 
  //if(LMAX>RMAX){
  //  grid.ngrid = (int)(cp.lbox/LMAX);
  //}else{ 
  //  grid.ngrid = (int)(cp.lbox/RMAX);
  //}

  grid.ngrid = (int)(cp.lbox/2000.);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }

  grid_init();
  grid_build();
  //////////////////////////////////////////////////

  c = (int *) malloc(NTHREADS*sizeof(int));
  pfdens = (FILE **) malloc(NTHREADS*sizeof(FILE));
  #ifdef CALCULA_MEDIA
    pfvel  = (FILE **) malloc(NTHREADS*sizeof(FILE));
  #endif
  #ifdef SAVEPART
    pfpart = (FILE **) malloc(NTHREADS*sizeof(FILE));
  #endif
  #ifdef EXTEND
    pfextend  = (FILE **) malloc(NTHREADS*sizeof(FILE));
  #endif
 
  for(i=0;i<NTHREADS;i++)
  {
    c[i] = 0;

    set_name(name,NNN,"densidad_pares");
    sprintf(filename,"%.2d_%s_%.2f_%.2f.%.2d.bin",snap.num,name,fof[0],fof[1],i);
    pfdens[i]=fopen(filename,"w");
    fwrite(&c[i],sizeof(int),1,pfdens[i]);        
    
    fprintf(stdout,"%s\n",filename);

    #ifdef CALCULA_MEDIA
      set_name(name,NNN,"vmedia_pares");
      sprintf(filename,"../%.2d_%s_%.2f_%.2f.%.2d.bin",snap.num,name,fof[0],fof[1],i);
      pfvel[i] = fopen(filename,"w");
      fwrite(&c[i],sizeof(int),1,pfvel[i]);        
      fwrite(&ncil,sizeof(int),1,pfvel[i]);
    #endif

    #ifdef SAVEPART
      set_name(name,NNN,"particulas_pares");
      sprintf(filename,"%.2d_%s_%.2f_%.2f.%.2d.bin",snap.num,name,fof[0],fof[1],i);
      pfpart[i] = fopen(filename,"w");
      fwrite(&c[i],sizeof(int),1,pfpart[i]);    
    #endif

    #ifdef EXTEND
      set_name(name,NNN,"extend_pares");
      sprintf(filename,"%.2d_%s_%.2f_%.2f.%.2d.bin",snap.num,name,fof[0],fof[1],i);
      pfextend[i]=fopen(filename,"w");
      fwrite(&c[i],sizeof(int),1,pfextend[i]);        
    #endif
  }

  #ifdef CALCULA_MEDIA

  BLUE("**** CALCULA VEL MEDIAS ******\n");
  #ifdef FIXED_SEPARATION

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(i,Tid, \
    rcil2,r,rsep,rlong,rlong_2) \
    shared(Gr,Seg,cp,ncil,binsper,c,pfvel,RLEN,RMAX,RMIN,nbins)

  #else

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(i,Tid, \
    rcil2,r,rsep,rlong,rlong_2) \
    shared(Gr,Seg,cp,ncil,binsper,c,pfvel,RMAX,RMIN,nbins)

  #endif
  for(i=0;i<cp.nseg;i++)
  {
    Tid = omp_get_thread_num();
    Seg[i].Vmedia = (type_real *) calloc(3*ncil,sizeof(type_real));
    rcil2  = (type_real *) malloc(binsper*sizeof(type_real));

    rsep    = Seg[i].len/(type_real)(nbins-1);
    rlong   = 2.0*rsep; 
    rlong_2 = rlong*0.5;

    #ifdef FIXED_SEPARATION
      r = RLEN; // RLEN<RMAX
    #else
      r = Seg[i].len;
    #endif

    if(r>RMAX)
    {
    #ifdef BIN_LOG
      logspace(rcil2,RMAX,RMIN,binsper);
    #else
      linspace(rcil2,RMAX,0.0,binsper);
    #endif
    }else{
    #ifdef BIN_LOG
      logspace(rcil2,r,r/50.,binsper);
    #else
      linspace(rcil2,r,0.0,binsper);
    #endif
    }

    cylinder_mean(i,rsep,rlong,rlong_2,rcil2);

    fwrite(&i,sizeof(int),1,pfvel[Tid]);
    fwrite(&Seg[i].Vmedia[0],sizeof(type_real),3*ncil,pfvel[Tid]);
    c[Tid]++;    
    free(rcil2);
  } //FINALIZA EL PARALELO
  BLUE("******************************\n");

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfvel[i]);
    fwrite(&c[i],sizeof(int),1,pfvel[i]);
    fclose(pfvel[i]);
    c[i] = 0;
  }

  #endif

  BLUE("********** CALCULA ***********\n");
  #ifdef SAVEPART
  
    #ifdef EXTEND 

      #ifdef FIXED_SEPARATION

        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfextend,pfdens,pfpart,Gr,P,Seg,cp, \
        RLEN,RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

      #else

        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfextend,pfdens,pfpart,Gr,P,Seg,cp, \
        RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

      #endif

    #else

      #ifdef FIXED_SEPARATION

        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfdens,pfpart,Gr,P,Seg,cp, \
        RLEN,RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)
    
      #else
  
        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfdens,pfpart,Gr,P,Seg,cp, \
        RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

      #endif

    #endif

  #else
  
    #ifdef EXTEND

      #ifdef FIXED_SEPARATION

        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfextend,pfdens,Gr,P,Seg,cp, \
        RLEN,RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

      #else

        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfextend,pfdens,Gr,P,Seg,cp, \
        RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

      #endif

    #else

      #ifdef FIXED_SEPARATION
      
        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfdens,Gr,P,Seg,cp, \
        RLEN,RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

      #else

        #pragma omp parallel for num_threads(NTHREADS) \
        schedule(dynamic) default(none) private(i,j,l,idim,rcil2,Tid, \
        r,rsep,rlong,rlong_2) shared(pfdens,Gr,P,Seg,cp, \
        RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)
    
      #endif

    #endif

  #endif
  for(i=0;i<cp.nseg;i++)
  {
    type_real  aux_mean, aux_rms;
    type_real *nodo;
    type_real *vmean_par, *vquad_par;
    type_real *vmean_perp, *vquad_perp;
    type_real **vmean, **vquad;
    struct node_sph **root;
    #ifdef SAVEPART
    int size_part = 0;
    int *vpart;
    vpart = (int *) calloc(cp.npart/32+1,sizeof(int));
    #endif

    rcil2  = (type_real *) malloc(binsper*sizeof(type_real));

    rsep    = Seg[i].len/(type_real)(nbins-1);
    rlong   = 2.0*rsep; 
    rlong_2 = rlong*0.5;

    #ifdef FIXED_SEPARATION
      r = RLEN; // RLEN<RMAX
    #else
      r = Seg[i].len;      
    #endif

    if(r>RMAX)
    {
      #ifdef BIN_LOG
        logspace(rcil2,RMAX,RMIN,binsper);
      #else
        linspace(rcil2,RMAX,0.0,binsper);
      #endif
    }else{
      #ifdef BIN_LOG
        logspace(rcil2,r,r/50.,binsper);
      #else
        linspace(rcil2,r,0.0,binsper);
      #endif
    }

    Tid = omp_get_thread_num();

    nodo = (type_real *) malloc(nbins*sizeof(type_real));
    vmean_par = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_par = (type_real *) malloc(ncil*sizeof(type_real));
    vmean_perp = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_perp = (type_real *) malloc(ncil*sizeof(type_real));   

    root = (struct node_sph **) malloc(ncil*sizeof(struct node_sph *));
    vmean = (type_real **) malloc(ncil*sizeof(type_real *));
    vquad = (type_real **) malloc(ncil*sizeof(type_real *));
    for(j=0;j<ncil;j++)
    {
      root[j] = (struct node_sph *) malloc(nbins*sizeof(struct node_sph));
      vmean[j] = (type_real *) malloc(3*sizeof(type_real));
      vquad[j] = (type_real *) malloc(3*sizeof(type_real));
    }

    #ifdef SAVEPART
      calcular_mean(i,binsper,rsep,rlong,rlong_2,rcil2,\
      vmean,vquad,vmean_par,vquad_par,\
      vmean_perp,vquad_perp,nodo,root,\
      &size_part,vpart);
    #else  
      calcular_mean(i,binsper,rsep,rlong,rlong_2,rcil2,\
      vmean,vquad,vmean_par,vquad_par,\
      vmean_perp,vquad_perp,nodo,root);
    #endif  

    fwrite(&i,sizeof(int),1,pfdens[Tid]);                   // Num Filamento
    fwrite(&Seg[i].flag,sizeof(int),1,pfdens[Tid]);         // escribe la bandera
    fwrite(&Seg[i].len,sizeof(type_real),1,pfdens[Tid]);    // escribe la longitud del segmento

    r = Gr[Seg[i].list[0]].NumPart*cp.Mpart;
    fwrite(&r,sizeof(type_real),1,pfdens[Tid]);             // escribe la masa del primer nodo
    r = Gr[Seg[i].list[Seg[i].size-1]].NumPart*cp.Mpart;
    fwrite(&r,sizeof(type_real),1,pfdens[Tid]);             // escribe la masa del ultimo nodo

    fwrite(&nbins,sizeof(int),1,pfdens[Tid]);               // la cantidad de cilindros que hizo

    for(j=0;j<nbins;j++)
    {
      //r = nodo[j]/Seg[i].len; // Para normalizarlo
      r = nodo[j];              // Sin normalizacion
      fwrite(&r,sizeof(type_real),1,pfdens[Tid]);
    }

    for(j=0;j<ncil;j++)
    {

      //r  = sqrt(rcil2[j]/rcil2[0]);
      //r += sqrt(rcil2[j+1]/rcil2[0]);
      //r *= 0.5;

      r  = sqrt(rcil2[j]);
      r += sqrt(rcil2[j+1]);
      r *= 0.5;

      fwrite(&r,sizeof(type_real),1,pfdens[Tid]);       // rbin_centre

      for(l=0;l<nbins;l++)
      {
        root[j][l].b -= 1.0;

        fwrite(&root[j][l].a,sizeof(int),1,pfdens[Tid]);       // npart
        fwrite(&root[j][l].b,sizeof(type_real),1,pfdens[Tid]); // pho - overdense

        ////////////////////////////////////////////////////////////////////////////////////

        for(idim=0;idim<3;idim++)
        {
          root[j][l].c[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
          root[j][l].d[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        }
        root[j][l].e *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].f *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].g *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].h *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;

        aux_mean = aux_rms = 0.0;

        for(idim=0;idim<3;idim++)
        {
          aux_mean  += root[j][l].c[idim];
          aux_rms   += root[j][l].d[idim] - root[j][l].c[idim]*root[j][l].c[idim];
        }
        
        aux_mean /= 3.0;        
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms))/3.0 : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean
        fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms

        aux_mean = root[j][l].e;
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].f - root[j][l].e*root[j][l].e)) : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean par
        fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms  par

        aux_mean = root[j][l].g;
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].h - root[j][l].g*root[j][l].g)) : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean perp        
        fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms  perp
      }

      free(vmean[j]);
      free(vquad[j]);
      free(root[j]);
    }

    free(nodo);
    free(root);
    free(vmean);
    free(vquad);
    free(vmean_par);
    free(vquad_par);
    free(vmean_perp);
    free(vquad_perp);

    #ifdef EXTEND

      #ifdef SAVEPART
        ext_mean(pfextend[Tid],i,binsper,rsep,rlong,rlong_2,rcil2,&size_part,vpart);
      #else
        ext_mean(pfextend[Tid],i,binsper,rsep,rlong,rlong_2,rcil2);
      #endif

    #endif

    #ifdef SAVEPART
    fwrite(&i,sizeof(int),1,pfpart[Tid]);
    fwrite(&size_part,sizeof(int),1,pfpart[Tid]);

    for(l=0;l<cp.npart;l++)
    {
      if(TestBit(vpart,l))    
        fwrite(&P[l].id,sizeof(int),1,pfpart[Tid]);
    }

    free(vpart);
    #endif

    c[Tid]++;
    totfil++;
    free(rcil2);

  } //FINALIZA EL PARALELO

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfdens[i]);
    fwrite(&c[i],sizeof(int),1,pfdens[i]);
    fclose(pfdens[i]);

    #ifdef EXTEND
      rewind(pfextend[i]);
      fwrite(&c[i],sizeof(int),1,pfextend[i]);
      fclose(pfextend[i]);
    #endif

    #ifdef SAVEPART
      rewind(pfpart[i]);
      fwrite(&c[i],sizeof(int),1,pfpart[i]);
      fclose(pfpart[i]);
    #endif
  }

  free(c);
  grid_free();

  BLUE("******************************\n");
  fprintf(stdout,"Calcula %d Filamentos\n",totfil);  
  BLUE("******************************\n");

  return;
}

#ifdef SAVEPART
  void calc_part(type_real * const Pos_cent, type_real * const Vmedia, int * vec, int * ncont, int * npart, type_real ** mean, type_real ** quad, \
  type_real * versor, type_real * mean_par, type_real * quad_par, type_real * mean_perp, type_real * quad_perp, type_real * const rcil2, type_real rlong_2, int binsper)
#else        
  void calc_part(type_real * const Pos_cent, type_real * const Vmedia, int * npart, type_real ** mean, type_real ** quad, \
  type_real * versor, type_real * mean_par, type_real * quad_par, type_real * mean_perp, type_real * quad_perp, type_real * const rcil2, type_real rlong_2, int binsper)
#endif
{
  int i,j,k;
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  type_real lbox,fac,lbox2;
  type_real dot,dis;
  type_real aux_real, aux_dot, aux_mod;
  int ngrid, ilim;
  type_real xxx[3];
  type_real Posprima[3],Vprima[3];
  #ifdef HUECOS
  int bin;
  #endif

  ////////////////////////////////////////////////////
  memset(npart,0,binsper*sizeof(int));
  memset(mean_par,0.0,binsper*sizeof(type_real));
  memset(quad_par,0.0,binsper*sizeof(type_real));
  memset(mean_perp,0.0,binsper*sizeof(type_real));
  memset(quad_perp,0.0,binsper*sizeof(type_real));
  for(j=0;j<binsper;j++)
  {
    memset(mean[j],0.0,3*sizeof(type_real));
    memset(quad[j],0.0,3*sizeof(type_real));
  }
  ////////////////////////////////////////////////////

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  dis = sqrt(rcil2[0]);
  ilim = dis>rlong_2 ? (int)(dis*fac)+1 : (int)(rlong_2*fac)+1;

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

          Posprima[0] = P[i].Pos[0] - Pos_cent[0];
          Posprima[1] = P[i].Pos[1] - Pos_cent[1];
          Posprima[2] = P[i].Pos[2] - Pos_cent[2];

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

            #ifdef HUECOS
            //if(dis<rcil2[0] && dis>rcil2[binsper])
            if(dis<rcil2[0])
            {

              bin = 0;
              for(j=1;j<binsper;j++)
              {
                bin = dis<rcil2[j] ? j : bin;
              }

              #ifdef SAVEPART
              if(bin>=3 && !TestBit(vec,i))
              {
                SetBit(vec,i);
                *ncont = *ncont + 1;
              }
              #endif

              aux_mod = aux_dot = 0.0;
              for(k=0;k<3;k++)
              {
                Vprima[k] = P[i].Vel[k] - Vmedia[3*bin+k];

                mean[bin][k] += Vprima[k];
                quad[bin][k] += (Vprima[k]*Vprima[k]);
  
                aux_dot += (Vprima[k]*versor[k]);
  
                xxx[k] = (Posprima[k] - dot*versor[k]);                               
  
                #ifdef PERIODIC
                if(xxx[k]>  lbox2) xxx[k] -= lbox;
                if(xxx[k]< -lbox2) xxx[k] += lbox;
                #endif
  
                aux_mod += (xxx[k]*xxx[k]);
              }
  
              mean_par[bin] += aux_dot;
              quad_par[bin] += (aux_dot*aux_dot);
              aux_mod = 1.0/sqrt(aux_mod);
  
              aux_real = 0.0;
              for(k=0;k<3;k++)
              {
                xxx[k] *= aux_mod;
                aux_real += (Vprima[k]*xxx[k]);
              }              
  
              mean_perp[bin] += aux_real;
              quad_perp[bin] += (aux_real*aux_real);
   
              npart[bin]++;
  
            }

            #else

            for(j=0;j<binsper;j++)
            {
              if(dis<rcil2[j])
              {
                if(j==0)
                {

                  #ifdef SAVEPART

                  if(!TestBit(vec,i))
                  {
                    SetBit(vec,i);
                    *ncont = *ncont + 1;
                  }
                  #endif

                  aux_mod = aux_dot = 0.0;
                  for(k=0;k<3;k++)
                  {
                    Vprima[k] = P[i].Vel[k] - Vmedia[j+binsper*k];

                    mean[j][k] += Vprima[k];
                    quad[j][k] += (Vprima[k]*Vprima[k]);

                    aux_dot += (Vprima[k]*versor[k]);

                    xxx[k] = (Posprima[k] - dot*versor[k]);                               

                    #ifdef PERIODIC
                    if(xxx[k]>  lbox2) xxx[k] -= lbox;
                    if(xxx[k]< -lbox2) xxx[k] += lbox;
                    #endif

                    aux_mod += (xxx[k]*xxx[k]);
                  }

                  mean_par[j] += aux_dot;
                  quad_par[j] += (aux_dot*aux_dot);
                  aux_mod = 1.0/sqrt(aux_mod);

                  aux_real = 0.0;
                  for(k=0;k<3;k++)
                  {
                    xxx[k] *= aux_mod;
                    aux_real += (Vprima[k]*xxx[k]);
                  }              

                  mean_perp[j] += aux_real;
                  quad_perp[j] += (aux_real*aux_real);
 
                }else{

                  for(k=0;k<3;k++)
                  {
                    mean[j][k] += Vprima[k];
                    quad[j][k] += (Vprima[k]*Vprima[k]);
                  }              

                  mean_par[j] += aux_dot;
                  quad_par[j] += (aux_dot*aux_dot);

                  mean_perp[j] += aux_real;
                  quad_perp[j] += (aux_real*aux_real); 
                  
                }

                npart[j]++;

              }else{

                break;

              }
            }

          #endif

          }// cierra el if

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

void set_name(char * name, const int NNN, char * prefix)
{

  #ifdef FIXED_SEPARATION
    sprintf(name,"%.4d_fixed",NNN);
  #else
    sprintf(name,"%.4d_varia",NNN);
  #endif

  #ifdef MCRITIC
    sprintf(name,"%s_%s_cut",name,prefix);
  #else
    sprintf(name,"%s_%s",name,prefix);
  #endif

  #ifdef BIN_LOG
    sprintf(name,"%s_LOG",name);
  #else
    sprintf(name,"%s_LIN",name);
  #endif

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

#ifdef SAVEPART

  void calcular_mean(const int i, const int binsper, const type_real rsep, \
           const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
  type_real **vmean, type_real **vquad, type_real *vmean_par, type_real *vquad_par,\
  type_real *vmean_perp, type_real *vquad_perp, type_real *nodo, struct node_sph **root,\
  int *size_part, int *vpart)

#else

  void calcular_mean(const int i, const int binsper, const type_real rsep, \
           const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
  type_real **vmean, type_real **vquad, type_real *vmean_par, type_real *vquad_par,\
  type_real *vmean_perp, type_real *vquad_perp, type_real *nodo, struct node_sph **root)

#endif
{
  int k,l,j,idim;
  int bin,lbin;
  int vsize[ncil];
  type_real vdir[3];
  type_real Pos_cent[3];
  type_real mod_v[3];
  type_real r,rbin,racum;
  type_real *volinv;                // En Kpc
  type_real dens, lenrbin;

  volinv = (type_real *) malloc(ncil*sizeof(type_real));

  r = 0.0f;
  for(j=binsper-2;j>=0;j--) // arranco en binsper-2 porque el primer cilindro esta completo
  {
    volinv[j] = M_PI*rlong*rcil2[j];

    #ifdef HUECOS
     volinv[j] -= r;

     r = M_PI*rlong*rcil2[j];
    #endif  

    volinv[j] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[j])/(type_real)cp.npart;
  }

  racum   = 0.0;
  lenrbin = 0.0;
  l = 0;

  ///////////////////////////////////////////////////////
  r = 0.0;
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
  ///////////////////////////////////////////////////////

  #ifdef SAVEPART
  calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
  #else
  calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
  #endif

  nodo[l] = lenrbin; 

  for(j=0;j<ncil;j++)
  {
    dens = (type_real)vsize[j]*volinv[j];
    push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
    vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
  }
  l++; // sumo

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

    #ifdef SAVEPART
    calc_part(Pos_cent,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #else
    calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #endif
    lenrbin += rsep;

    nodo[l] = lenrbin; 
    for(j=0;j<ncil;j++)
    {
      dens = (type_real)vsize[j]*volinv[j];
      push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
      vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
    }
    l++; // sumo

    r = 0.0;
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

      #ifdef SAVEPART
      calc_part(Pos_cent,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #else
      calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #endif

      lenrbin += rsep;
      nodo[l] = lenrbin; 

      for(j=0;j<ncil;j++)
      {
        dens = (type_real)vsize[j]*volinv[j];
        push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
        vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
      }

      l++; // sumo
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
  ///////////////////////////////////////////////////////

  #ifdef SAVEPART
  calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vpart,size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
  #else
  calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
  #endif

  nodo[l] = Seg[i].len; 
  for(j=0;j<ncil;j++)
  {
    dens = (type_real)vsize[j]*volinv[j];
    push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
    vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
  }
  l++; // sumo

  assert(l==nbins);

  free(volinv);
}

#ifdef CALCULA_MEDIA

void cylinder_mean(const int i, const type_real rsep, const type_real rlong, \
                   const type_real rlong_2, type_real * const rcil2)
{
  int *npart;
  int k,l,idim;
  int bin,lbin;
  type_real vdir[3];
  type_real Pos_cent[3];
  type_real mod_v[3];
  type_real r,rbin,racum;

  npart = (int *) calloc(ncil,sizeof(int));
  l = 0;
  racum = 0.0f;

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
    vdir[idim] *= (1.0/r);
  ///////////////////////////////////////////////////////

  calc_media(Gr[Seg[i].list[0]].Pos,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);             
  l++;

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

    calc_media(Pos_cent,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);             
    l++;

    r = 0.0;
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
      calc_media(Pos_cent,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);

      l++;
    }

    racum = r-(type_real)lbin*rsep;
  }

  if(l==nbins)
    l--;     // sobre escribo el ultimo

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
  ///////////////////////////////////////////////////////

  calc_media(Gr[Seg[i].list[Seg[i].size-1]].Pos,vdir,npart,Seg[i].Vmedia,rcil2,rlong_2,ncil);             
  l++;

  assert(l==nbins);

  for(bin=0;bin<ncil;bin++)
  {
    Seg[i].Vmedia[3*bin + 0] *= (npart[bin]>0 ? 1.0f/(type_real)npart[bin] : 0.0);
    Seg[i].Vmedia[3*bin + 1] *= (npart[bin]>0 ? 1.0f/(type_real)npart[bin] : 0.0);
    Seg[i].Vmedia[3*bin + 2] *= (npart[bin]>0 ? 1.0f/(type_real)npart[bin] : 0.0);
  }

  free(npart);

}

void calc_media(type_real * const Pos_cent, type_real * const versor, int *numpart, type_real *vel_media, type_real * const rcil2, const type_real rlong_2,
const int binsper)
{
  int i,j,bin;
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  type_real Posprima[3];
  type_real lbox,fac,lbox2;
  type_real dot,dis;
  int ngrid, ilim;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  dis = sqrt(rcil2[0]);
  ilim = dis>rlong_2 ? (int)(dis*fac)+1 : (int)(rlong_2*fac)+1;

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
          Posprima[0] = P[i].Pos[0] - Pos_cent[0];
          Posprima[1] = P[i].Pos[1] - Pos_cent[1];
          Posprima[2] = P[i].Pos[2] - Pos_cent[2];

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

            //if(dis<rcil2[0] && dis>rcil2[binsper])
            if(dis<rcil2[0])
            {
              bin = 0;
              for(j=1;j<binsper;j++)
              {
                bin = dis<rcil2[j] ? j : bin;
              }

              vel_media[3*bin+0] += P[i].Vel[0];
              vel_media[3*bin+1] += P[i].Vel[1];
              vel_media[3*bin+2] += P[i].Vel[2];
              numpart[bin]++;

            }
          } // cierra dis

          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

#endif

#ifdef EXTEND

  #ifdef SAVEPART
  
    void ext_mean(FILE *pfextend; const int i, const int binsper, const type_real rsep,\
    const type_real rlong, const type_real rlong_2, type_real * const rcil2,\
    int *size_part, int *vpart)
 
  #else
  
    void ext_mean(FILE *pfextend, const int i, const int binsper, const type_real rsep,\
    const type_real rlong, const type_real rlong_2, type_real * const rcil2)
 
  #endif
  {
    int l,j,idim,lbin;
    int vsize[ncil];
    int *npart;
    type_real aux_mean, aux_rms;
    type_real r, dens;
    type_real vdir[3];
    type_real Pos_cent[3];
    type_real *volinv;                // En Kpc
    //type_real *Vmean;
    type_real *nodo;
    type_real *vmean_par, *vquad_par;
    type_real *vmean_perp, *vquad_perp;
    type_real **vmean, **vquad;
    struct node_sph **root;
  
    ////////////////////////////////////////////////////////////////////////////////
    l = 0;
    lbin = nbins%2==0 ? nbins/2 : (int)floor((float)nbins/2);
    assert(lbin!=0);

    volinv = (type_real *) malloc(ncil*sizeof(type_real));
    //Vmean  = (type_real *) calloc(3*ncil,sizeof(type_real));
    npart  = (int *) calloc(ncil,sizeof(int));

    nodo = (type_real *) malloc(lbin*sizeof(type_real));
    vmean_par = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_par = (type_real *) malloc(ncil*sizeof(type_real));
    vmean_perp = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_perp = (type_real *) malloc(ncil*sizeof(type_real));   

    root = (struct node_sph **) malloc(ncil*sizeof(struct node_sph *));
    vmean = (type_real **) malloc(ncil*sizeof(type_real *));
    vquad = (type_real **) malloc(ncil*sizeof(type_real *));
    for(j=0;j<ncil;j++)
    {
      root[j] = (struct node_sph *) malloc(lbin*sizeof(struct node_sph));
      vmean[j] = (type_real *) malloc(3*sizeof(type_real));
      vquad[j] = (type_real *) malloc(3*sizeof(type_real));
    }
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    j = 2*lbin;
    fwrite(&i,sizeof(int),1,pfextend);                   // Num Filamento
    fwrite(&Seg[i].len,sizeof(type_real),1,pfextend);    // escribe la longitud del segmento
    fwrite(&j,sizeof(int),1,pfextend);                   // la cantidad de cilindros que voy a hacer
    ////////////////////////////////////////////////////////////////////////////////      

    ////////////////////////////////////////////////////////////////////////////////
    r = 0.0f;
    for(j=binsper-2;j>=0;j--) // arranco en binsper-2 porque el primer cilindro esta completo
    {
      volinv[j] = M_PI*rlong*rcil2[j];
  
      #ifdef HUECOS
       volinv[j] -= r;
  
       r = M_PI*rlong*rcil2[j];
      #endif  
  
      volinv[j] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[j])/(type_real)cp.npart;
    }
    ////////////////////////////////////////////////////////////////////////////////

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
      vdir[idim] *= (1/r); // directo

    /*
    #ifdef CALCULA_MEDIA

    ////////////////////////////////////////////////////////////////////////////////

    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] = Gr[Seg[i].list[0]].Pos[idim] - (type_real)(lbin-1)*rsep*vdir[idim];
      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }
    ///////////////////////////////////////////////////////////////////////////////
    calc_media(Pos_cent,vdir,npart,Vmean,rcil2,rlong_2,ncil);
    l++;

    while(l<lbin)
    {
      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }
      
      if(l==lbin-1)
        calc_media(Gr[Seg[i].list[0]].Pos,vdir,npart,Vmean,rcil2,rlong_2,ncil);
      else
        calc_media(Pos_cent,vdir,npart,Vmean,rcil2,rlong_2,ncil);
      l++;
    }    
 
    assert(l==lbin);

    for(l=0;l<ncil;l++)
    {
      Vmean[l]        *= (npart[l]>0 ? 1.0f/(type_real)npart[l] : 0.0);
      Vmean[l+ncil]   *= (npart[l]>0 ? 1.0f/(type_real)npart[l] : 0.0);
      Vmean[l+2*ncil] *= (npart[l]>0 ? 1.0f/(type_real)npart[l] : 0.0);
    }
    ////////////////////////////////////////////////////////////////////////////////

    #endif
    */

    l = 0;
    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] = Gr[Seg[i].list[0]].Pos[idim] - (lbin-1)*rsep*vdir[idim];
      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }

    #ifdef SAVEPART
      //calc_part(Pos_cent,Vmean,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
      calc_part(Pos_cent,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
    #else
      //calc_part(Pos_cent,Vmean,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #endif

    nodo[l] = -(type_real)(lbin-1)*rsep;
    for(j=0;j<ncil;j++)
    {
      dens = (type_real)vsize[j]*volinv[j];
      push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
      vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
    }
    l++; // sumo  
 
    while(l<lbin)
    {  
      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }

      #ifdef SAVEPART
        if(l==lbin-1)
          //calc_part(Gr[Seg[i].list[0]].Pos,Vmean,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
          calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
        else
          //calc_part(Pos_cent,Vmean,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
          calc_part(Pos_cent,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #else
        if(l==lbin-1)
          //calc_part(Gr[Seg[i].list[0]].Pos,Vmean,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
          calc_part(Gr[Seg[i].list[0]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
        else
          //calc_part(Pos_cent,Vmean,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
          calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #endif
 
      nodo[l] = nodo[l-1]+rsep;   
      for(j=0;j<ncil;j++)
      {
        dens = (type_real)vsize[j]*volinv[j];
        push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
        vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
      }
  
      l++; // sumo
    }
  
    assert(l==lbin);
    
    ////////////////////////////////////////////////////////////////////////////////
    for(j=0;j<lbin;j++)
    {
      nodo[j] = j==lbin-1 ? 0.0f : nodo[j];
      //r = nodo[j]/Seg[i].len; // Para normalizarlo
      r = nodo[j];              // Sin normalizacion
      fwrite(&r,sizeof(type_real),1,pfextend);
    }

    for(j=0;j<ncil;j++)
    {

      r  = sqrt(rcil2[j]);
      r += sqrt(rcil2[j+1]);
      r *= 0.5;

      fwrite(&r,sizeof(type_real),1,pfextend);       // rbin_centre

      for(l=0;l<lbin;l++)
      {
        root[j][l].b -= 1.0;

        fwrite(&root[j][l].a,sizeof(int),1,pfextend);       // npart
        fwrite(&root[j][l].b,sizeof(type_real),1,pfextend); // pho - overdense

        ////////////////////////////////////////////////////////////////////////////////////

        for(idim=0;idim<3;idim++)
        {
          root[j][l].c[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
          root[j][l].d[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        }
        root[j][l].e *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].f *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].g *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].h *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;

        aux_mean = aux_rms = 0.0;

        for(idim=0;idim<3;idim++)
        {
          aux_mean  += root[j][l].c[idim];
          aux_rms   += root[j][l].d[idim] - root[j][l].c[idim]*root[j][l].c[idim];
        }
        
        aux_mean /= 3.0;        
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms))/3.0 : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfextend); // mean
        fwrite(&aux_rms,sizeof(type_real),1,pfextend);  // rms

        aux_mean = root[j][l].e;
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].f - root[j][l].e*root[j][l].e)) : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfextend); // mean par
        fwrite(&aux_rms,sizeof(type_real),1,pfextend);  // rms  par

        aux_mean = root[j][l].g;
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].h - root[j][l].g*root[j][l].g)) : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfextend); // mean perp        
        fwrite(&aux_rms,sizeof(type_real),1,pfextend);  // rms  perp
      }

      free(vmean[j]);
      free(vquad[j]);
      free(root[j]);
    }

    free(nodo);
    free(root);
    free(vmean);
    free(vquad);
    free(vmean_par);
    free(vquad_par);
    free(vmean_perp);
    free(vquad_perp);
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    l = 0;
    //memset(Vmean,0.0,3*ncil*sizeof(type_real));
    memset(npart,0,ncil*sizeof(int));
    nodo = (type_real *) malloc(lbin*sizeof(type_real));
    vmean_par = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_par = (type_real *) malloc(ncil*sizeof(type_real));
    vmean_perp = (type_real *) malloc(ncil*sizeof(type_real));
    vquad_perp = (type_real *) malloc(ncil*sizeof(type_real));   

    root = (struct node_sph **) malloc(ncil*sizeof(struct node_sph *));
    vmean = (type_real **) malloc(ncil*sizeof(type_real *));
    vquad = (type_real **) malloc(ncil*sizeof(type_real *));
    for(j=0;j<ncil;j++)
    {
      root[j] = (struct node_sph *) malloc(lbin*sizeof(struct node_sph));
      vmean[j] = (type_real *) malloc(3*sizeof(type_real));
      vquad[j] = (type_real *) malloc(3*sizeof(type_real));
    }
    ////////////////////////////////////////////////////////////////////////////////

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
      vdir[idim] *= (1/r); // reversa

    /*
    #ifdef CALCULA_MEDIA

    ////////////////////////////////////////////////////////////////////////////////
    for(idim=0;idim<3;idim++)
    {
      Pos_cent[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim];
      #ifdef PERIODIC
      Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
      Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
      #endif
    }
    ///////////////////////////////////////////////////////////////////////////////

    calc_media(Gr[Seg[i].list[Seg[i].size-1]].Pos,vdir,npart,Vmean,rcil2,rlong_2,ncil);             
    l++;
   
    while(l<lbin)
    {
      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }
 
      calc_media(Pos_cent,vdir,npart,Vmean,rcil2,rlong_2,ncil);
      l++;
    }

    assert(l==lbin);

    for(l=0;l<ncil;l++)
    {
      Vmean[l]        *= (npart[l]>0 ? 1.0f/(type_real)npart[l] : 0.0);
      Vmean[l+ncil]   *= (npart[l]>0 ? 1.0f/(type_real)npart[l] : 0.0);
      Vmean[l+2*ncil] *= (npart[l]>0 ? 1.0f/(type_real)npart[l] : 0.0);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    #endif
    */

    l = 0;
    for(idim=0;idim<3;idim++)
      Pos_cent[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim];

    #ifdef SAVEPART
    //calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Vmean,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
    calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);
    #else
    //calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Vmean,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    calc_part(Gr[Seg[i].list[Seg[i].size-1]].Pos,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
    #endif
  
    nodo[l] = 0.0f;   
    for(j=0;j<ncil;j++)
    {
      dens = (type_real)vsize[j]*volinv[j];
      push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
      vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
    }
    l++; // sumo  
 
    while(l<lbin)
    {  
      for(idim=0;idim<3;idim++)
      {
        Pos_cent[idim] += rsep*vdir[idim];
        #ifdef PERIODIC
        Pos_cent[idim] = Pos_cent[idim] >= cp.lbox ? Pos_cent[idim]-cp.lbox : Pos_cent[idim];
        Pos_cent[idim] = Pos_cent[idim] <      0.0 ? Pos_cent[idim]+cp.lbox : Pos_cent[idim];
        #endif
      }
  
      #ifdef SAVEPART
        //calc_part(Pos_cent,Vmean,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
        calc_part(Pos_cent,Seg[i].Vmedia,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #else
        //calc_part(Pos_cent,Vmean,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
        calc_part(Pos_cent,Seg[i].Vmedia,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #endif
  
      nodo[l] = nodo[l-1]+rsep;   
      for(j=0;j<ncil;j++)
      {
        dens = (type_real)vsize[j]*volinv[j];
        push_array(&root[j][l],vsize[j],dens,vmean[j],vquad[j], \
        vmean_par[j],vquad_par[j],vmean_perp[j],vquad_perp[j]);
      }
  
      l++; // sumo
    }
  
    assert(l==lbin);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(j=0;j<lbin;j++)
    {
      nodo[j] = j==0 ? Seg[i].len : Seg[i].len+nodo[j];
      //r = nodo[j]/Seg[i].len; // Para normalizarlo
      r = nodo[j];              // Sin normalizacion
      fwrite(&r,sizeof(type_real),1,pfextend);
    }

    for(j=0;j<ncil;j++)
    {

      r  = sqrt(rcil2[j]);
      r += sqrt(rcil2[j+1]);
      r *= 0.5;

      fwrite(&r,sizeof(type_real),1,pfextend);       // rbin_centre

      for(l=0;l<lbin;l++)
      {
        root[j][l].b -= 1.0;

        fwrite(&root[j][l].a,sizeof(int),1,pfextend);       // npart
        fwrite(&root[j][l].b,sizeof(type_real),1,pfextend); // pho - overdense

        ////////////////////////////////////////////////////////////////////////////////////

        for(idim=0;idim<3;idim++)
        {
          root[j][l].c[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
          root[j][l].d[idim] *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        }
        root[j][l].e *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].f *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].g *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;
        root[j][l].h *= (root[j][l].a>0) ? 1./(type_real)root[j][l].a : 0.0;

        aux_mean = aux_rms = 0.0;

        for(idim=0;idim<3;idim++)
        {
          aux_mean  += root[j][l].c[idim];
          aux_rms   += root[j][l].d[idim] - root[j][l].c[idim]*root[j][l].c[idim];
        }
        
        aux_mean /= 3.0;        
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(aux_rms))/3.0 : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfextend); // mean
        fwrite(&aux_rms,sizeof(type_real),1,pfextend);  // rms

        aux_mean = root[j][l].e;
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].f - root[j][l].e*root[j][l].e)) : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfextend); // mean par
        fwrite(&aux_rms,sizeof(type_real),1,pfextend);  // rms  par

        aux_mean = root[j][l].g;
        aux_rms  = root[j][l].a>10 ? sqrt(fabs(root[j][l].h - root[j][l].g*root[j][l].g)) : 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfextend); // mean perp        
        fwrite(&aux_rms,sizeof(type_real),1,pfextend);  // rms  perp
      }

      free(vmean[j]);
      free(vquad[j]);
      free(root[j]);
    }

    free(nodo);
    free(root);
    free(vmean);
    free(vquad);
    free(vmean_par);
    free(vquad_par);
    free(vmean_perp);
    free(vquad_perp);
    //free(Vmean);  
    free(volinv);
    free(npart);
  }

#endif

#ifdef VEL_RELATIVA

void relativa_write(int NNN, type_real *fof)
{
  int i, idim;
  type_real vdir[3];
  type_real r;
  FILE *pfrela;
  char filename[200], name[200];

  BLUE("******************************\n");

  sprintf(filename,"Escribe el archivo para calcular las velocidades relativas\n");GREEN(filename);

  set_name(name,NNN,"vrelativa_pares");
  sprintf(filename,"../%.2d_%s_%.2f_%.2f.bin",snap.num,name,fof[0],fof[1]);
  pfrela = fopen(filename,"w");
  fwrite(&cp.nseg,sizeof(int),1,pfrela);        
 
  ///////////////////////////////////////////////////////
  
  for(i=0;i<cp.nseg;i++)
  {
    r = 0.0f;
    for(idim=0;idim<3;idim++)
    {
      vdir[idim] = Gr[Seg[i].list[Seg[i].size-1]].Pos[idim]-Gr[Seg[i].list[0]].Pos[idim];    
      #ifdef PERIODIC
      if(vdir[idim]> 0.5*cp.lbox) vdir[idim] -= cp.lbox;
      if(vdir[idim]<-0.5*cp.lbox) vdir[idim] += cp.lbox;
      #endif
      r += vdir[idim]*vdir[idim];
    }
    r = sqrt(r);

    for(idim=0;idim<3;idim++)
      vdir[idim] *= (1.0/r);

    fwrite(&i,sizeof(int),1,pfrela);
    fwrite(&Seg[i].flag,sizeof(int),1,pfrela);    
    fwrite(&Gr[Seg[i].list[0]].Vnod[0],sizeof(type_real),3,pfrela);
    fwrite(&Gr[Seg[i].list[Seg[i].size-1]].Vnod[0],sizeof(type_real),3,pfrela);
    fwrite(&vdir[0],sizeof(type_real),3,pfrela);

  } //FINALIZA EL PARALELO

  fclose(pfrela);

  ///////////////////////////////////////////////////////

  BLUE("******************************\n");
}

#endif
