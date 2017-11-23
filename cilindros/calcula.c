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
  char filename[200];
  int i, j, k, l, idim, Tid, *c;
  int totfil = 0;
  int binsper = ncil + 1;
  int bin,lbin;
  type_real *rcil2;                 // En Kpc                    
  type_real *volinv;                // En Kpc
  type_real rsep;                   // En Kpc
  type_real rlong;                  // En Kpc
  type_real rlong_2;                // En Kpc
  type_real RMAX, RMIN, LMAX;
  type_real vdir[3];
  type_real Pos_cent[3];
  type_real mod_v[3];
  type_real r,rbin,racum;
  type_real dens,lenrbin;
  FILE **pfdens;
  #ifdef CALCULA_MEDIA
  int   *npart;
  FILE **pfvel;
  #endif
  #ifdef SAVEPART
  FILE **pfpart;
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

  fprintf(stdout,"RMAX %f Mpc\n",RMAX/1000.);
  fprintf(stdout,"RMIN %f Mpc\n",RMIN/1000.);
  fprintf(stdout,"LMAX %f Mpc\n",LMAX/1000.);

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

  for(i=0;i<NTHREADS;i++)
  {
    c[i] = 0;

    #ifdef MCRITIC
    sprintf(filename,"%.2d_%.4d_densidad_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NNN,m_critica,fof[0],fof[1],i);
    #else
    sprintf(filename,"%.2d_%.4d_densidad_%.2f_%.2f.%.2d.bin",snap.num,NNN,fof[0],fof[1],i);
    #endif
    pfdens[i]=fopen(filename,"w");
    fwrite(&c[i],sizeof(int),1,pfdens[i]);        

    #ifdef CALCULA_MEDIA
      #ifdef MCRITIC
        sprintf(filename,"../%.2d_%.4d_vmedia_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NNN,m_critica,fof[0],fof[1],i);
      #else
        sprintf(filename,"../%.2d_%.4d_vmedia_%.2f_%.2f.%.2d.bin",snap.num,NNN,fof[0],fof[1],i);
      #endif
    pfvel[i] = fopen(filename,"w");
    fwrite(&c[i],sizeof(int),1,pfvel[i]);        
    fwrite(&ncil,sizeof(int),1,pfvel[i]);
    #endif

    #ifdef SAVEPART
       #ifdef MCRITIC
         sprintf(filename,"%.2d_%.4d_particulas_cut_%.2f_%.2f_%.2f.%.2d.bin",snap.num,NNN,m_critica,fof[0],fof[1],i);
       #else
         sprintf(filename,"%.2d_%.4d_particulas_%.2f_%.2f.%.2d.bin",snap.num,NNN,fof[0],fof[1],i);
       #endif
       pfpart[i] = fopen(filename,"w");
       fwrite(&c[i],sizeof(int),1,pfpart[i]);    
    #endif

  }

  #ifdef CALCULA_MEDIA

  BLUE("**** CALCULA VEL MEDIAS ******\n");
  #pragma omp parallel for num_threads(NTHREADS) \
  schedule(dynamic) default(none) private(i,idim,k,l,Tid,      \
  rcil2,vdir,Pos_cent,mod_v,r,rbin,racum,bin,lbin,npart,rsep,  \
  rlong,rlong_2) shared(Gr,Seg,cp,ncil,binsper,c,pfvel,RMAX,RMIN,nbins)
  for(i=0;i<cp.nseg;i++)
  {
    Tid = omp_get_thread_num();
    npart = (int *) calloc(ncil,sizeof(int));
    Seg[i].Vmedia = (type_real *) calloc(3*ncil,sizeof(type_real));
    rcil2  = (type_real *) malloc(binsper*sizeof(type_real));
    l = 0;

    rsep    = Seg[i].len/(type_real)nbins;      
    rlong   = 2.0*rsep; 
    rlong_2 = rlong*0.5;

    r = 0.5*Seg[i].len;

    if(r>RMAX)
    {
      logspace(rcil2,RMAX,RMIN,binsper);
    }else{
      logspace(rcil2,r,r/50.,binsper);
    }

    racum = 0.0;

    for(k=1;k<Seg[i].size;k++)
    {

      r = 0.0;
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
        vdir[idim] /= r;
        Pos_cent[idim] = Gr[Seg[i].list[k-1]].Pos[idim];
      }

      if(k==1)
      {
        rbin = 0.5*rsep;
      }else{
        rbin = rsep-racum;
      }

      if(r<rbin)
      {
        if(k==1) r += rbin;
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

    assert(l==nbins);

    for(bin=0;bin<ncil;bin++)
    {
      if(npart[bin]==0)
        continue;

      for(idim=0;idim<3;idim++)
        Seg[i].Vmedia[bin+ncil*idim] /= (type_real)npart[bin];
    }

    fwrite(&i,sizeof(int),1,pfvel[Tid]);
    fwrite(&Seg[i].Vmedia[0],sizeof(type_real),3*ncil,pfvel[Tid]);
    c[Tid]++;    
    free(npart);
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

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(i,j,k,l,idim,rcil2,Tid, \
    bin,lbin,volinv,vdir,Pos_cent,mod_v,r,rbin,lenrbin, \
    racum,dens,rsep,rlong,rlong_2) shared(pfdens,pfpart,Gr,P,Seg,cp, \
    RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

  #else

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(i,j,k,l,idim,rcil2,Tid, \
    bin,lbin,volinv,vdir,Pos_cent,mod_v,r,rbin,lenrbin, \
    racum,dens,rsep,rlong,rlong_2) shared(pfdens,Gr,P,Seg,cp, \
    RMAX,RMIN,c,binsper,ncil,nbins,stdout) reduction(+:totfil)

  #endif

  for(i=0;i<cp.nseg;i++)
  {
    int        vsize[binsper];
    type_real  aux_mean, aux_rms;
    type_real **vmean, **vquad;
    type_real *vmean_par, *vquad_par;
    type_real *vmean_perp, *vquad_perp;
    type_real *nodo;
    struct node_sph **root;
    #ifdef HUECOS
    type_real tmpvol;
    #endif
    #ifdef SAVEPART
    int   size_part = 0;
    int   *vpart;
    vpart = (int *) calloc(cp.npart/32+1,sizeof(int));
    #endif

    rcil2  = (type_real *) malloc(binsper*sizeof(type_real));
    volinv = (type_real *) malloc(binsper*sizeof(type_real));

    rsep = Seg[i].len/(type_real)nbins;      
    rlong = 2.0*rsep; 
    rlong_2= rlong*0.5;

    r = 0.5*Seg[i].len;

    if(r>RMAX)
    {
      //logspace(rcil2,RMAX,RMIN,binsper);
      logspace(rcil2,RMAX,RMIN,binsper);
    }else{
      logspace(rcil2,r,r/50.,binsper);
    }

    for(k=binsper-1;k>=0;k--)
    {
      #ifdef HUECOS
        volinv[k] = M_PI*rlong*rcil2[k];

        if(k!=binsper-1) 
          volinv[k] -= tmpvol;        

        tmpvol = M_PI*rlong*rcil2[k];
      #endif  
        volinv[k] = ((cp.lbox*cp.lbox*cp.lbox)/volinv[k])/(type_real)cp.npart;
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

    racum   = 0.0;
    lenrbin = 0.0;
    l = 0;

    for(k=1;k<Seg[i].size;k++)
    {
      
      r = 0.0;
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
        vdir[idim] /= r;
        Pos_cent[idim] = Gr[Seg[i].list[k-1]].Pos[idim];
      }

      if(k==1)
      {
        rbin = 0.5*rsep;
      }else{
        rbin = rsep-racum;
      }

      if(r<rbin)
      {
        if(k==1) r += rbin;
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
      calc_part(Pos_cent,i,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #else
      calc_part(Pos_cent,i,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
      #endif

      if(l==0)
      {       
        lenrbin += 0.5*rsep;
      }else{
        lenrbin += rsep;
      }

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
        calc_part(Pos_cent,i,vpart,&size_part,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
        #else
        calc_part(Pos_cent,i,vsize,vmean,vquad,vdir,vmean_par,vquad_par,vmean_perp,vquad_perp,rcil2,rlong_2,ncil);             
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

    assert(l==nbins);

    fwrite(&i,sizeof(int),1,pfdens[Tid]);                   // Num Filamento
    fwrite(&Seg[i].flag,sizeof(int),1,pfdens[Tid]);         // escribe la bandera
    fwrite(&Seg[i].len,sizeof(type_real),1,pfdens[Tid]);    // escribe la longitud del segmento

    r = Gr[Seg[i].list[0]].NumPart*cp.Mpart;
    fwrite(&r,sizeof(type_real),1,pfdens[Tid]);             // escribe la masa del primer nodo
    r = Gr[Seg[i].list[Seg[i].size-1]].NumPart*cp.Mpart;
    fwrite(&r,sizeof(type_real),1,pfdens[Tid]);             // escribe la masa del ultimo nodo

    fwrite(&nbins,sizeof(int),1,pfdens[Tid]);               // la cantidad de pelotas que hizo

    for(j=0;j<nbins;j++)
    {
      //r = nodo[j]/Seg[i].len; // Para normalizarlo
      r = nodo[j];              // Sin normalizacion
      fwrite(&r,sizeof(type_real),1,pfdens[Tid]);
    }

    for(j=0;j<ncil;j++)
    {

      r  = sqrt(rcil2[j]/rcil2[0]);
      r += sqrt(rcil2[j+1]/rcil2[0]);
      r *= 0.5;

      fwrite(&r,sizeof(type_real),1,pfdens[Tid]);       // rbin_centre

      for(l=0;l<nbins;l++)
      {

        root[j][l].b -= 1.0;

        fwrite(&root[j][l].a,sizeof(int),1,pfdens[Tid]);       // npart
        fwrite(&root[j][l].b,sizeof(type_real),1,pfdens[Tid]); // pho - overdense

        ////////////////////////////////////////////////////////////////////////////////////

        if(root[j][l].a!=0)
        {
          for(idim=0;idim<3;idim++)
          {
            root[j][l].c[idim] /= (type_real)root[j][l].a;
            root[j][l].d[idim] /= (type_real)root[j][l].a;
          }
          root[j][l].e /= (type_real)root[j][l].a;
          root[j][l].f /= (type_real)root[j][l].a;
          root[j][l].g /= (type_real)root[j][l].a;
          root[j][l].h /= (type_real)root[j][l].a;
        }

        aux_mean = aux_rms = 0.0;

        for(idim=0;idim<3;idim++)
        {
          aux_mean  += root[j][l].c[idim];
          aux_rms   += root[j][l].d[idim] - root[j][l].c[idim]*root[j][l].c[idim];
        }
        
        aux_mean /= 3.0;        
        if(root[j][l].a>10)
          aux_rms  = sqrt(fabs(aux_rms))/3.0;
        else
          aux_rms  = 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean
        fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms

        aux_mean = root[j][l].e;
        if(root[j][l].a>10)
          aux_rms  = sqrt(fabs(root[j][l].f - root[j][l].e*root[j][l].e));
        else
          aux_rms  = 0.0;

        fwrite(&aux_mean,sizeof(type_real),1,pfdens[Tid]); // mean par
        fwrite(&aux_rms,sizeof(type_real),1,pfdens[Tid]);  // rms  par

        aux_mean = root[j][l].g;
        if(root[j][l].a>10)
          aux_rms  = sqrt(fabs(root[j][l].h - root[j][l].g*root[j][l].g));
        else
         aux_rms  = 0.0;

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
    free(volinv);

  } //FINALIZA EL PARALELO

  for(i=0;i<NTHREADS;i++)
  {
    fprintf(stdout,"Tid %d Nfil %d\n",i,c[i]);  
    rewind(pfdens[i]);
    fwrite(&c[i],sizeof(int),1,pfdens[i]);
    fclose(pfdens[i]);

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
  void calc_part(type_real *Pos_cent, int icent, int *vec, int *ncont, int *npart, type_real **mean, type_real **quad, \
  type_real *versor, type_real *mean_par, type_real *quad_par, type_real *mean_perp, type_real *quad_perp, type_real *rcil2, type_real rlong_2, int binsper)
#else        
  void calc_part(type_real *Pos_cent, int icent, int *npart, type_real **mean, type_real **quad, \
  type_real *versor, type_real *mean_par, type_real *quad_par, type_real *mean_perp, type_real *quad_perp, type_real *rcil2, type_real rlong_2, int binsper)
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

          //if(point_inside(dot,rlong_2)==1)
          if(point_inside(dot,rlong_2) & 1)
          {

            dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2] - dot*dot;

            #ifdef HUECOS
            if(dis<rcil2[0] && dis>rcil2[binsper])
            //if(dis<rcil2[0])
            {

              bin = 0;
              for(j=1;j<binsper;j++)
              {
                if(dis<rcil2[j])
                 bin = j;
                else
                  break;
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
                Vprima[k] = P[i].Vel[k] - Seg[icent].Vmedia[bin+binsper*k];

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
                    Vprima[k] = P[i].Vel[k] - Seg[icent].Vmedia[j+binsper*k];

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

#ifdef CALCULA_MEDIA

void calc_media(type_real *Pos_cent, type_real *versor, int *numpart, type_real *vel_media, type_real *rcil2, type_real rlong_2, int binsper)
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

          //if(point_inside(dot,rlong_2)==1)
          if(point_inside(dot,rlong_2) & 1)
          {
            dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2] - dot*dot;

            //if(dis<rcil2[0] && dis>rcil2[binsper])
            if(dis<rcil2[0])
            {

             for(k=0;k<3;k++)
               vel_media[binsper*k] += P[i].Vel[k];

             numpart[0]++;

             for(j=1;j<binsper;j++)
              {
                if(dis<rcil2[j])
                {
                  for(k=0;k<3;k++)
                    vel_media[j+binsper*k] += P[i].Vel[k];

                  numpart[j]++;

                }else{

                  break;

                }
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

#endif

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
      // if there is an even number of elements, return mean of the two elements in the middle
      median = (a[j/2] + a[j/2 - 1])*0.5;
      fprintf(stdout,"median %f\n",median);
  } else {
      // else return the element in the middle
      median = a[j/2];
      fprintf(stdout,"median %f\n",median);
  }


  if(mean>median)
    *MAX = mean;
  else
    *MAX = median;

  //*MAX = 500.*ceil(*MAX/1000.);
  *MAX = 14000.;  
  *MIN = *MAX/50.;


  free(a);
  
  return;
}

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
