#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "iden.h"
#include "colores.h"
#include "peano.h"
#include "grid.h"

void Write_Groups(double fof);

int main(int argc, char **argv)
{
  int    i;
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  /*Lee archivos de la simulacion*/
  read_gadget();

  // Cambia origen de coordenadas y calcula la Peano Hilbert
  peano_hilbert();

  for(i=0;i<cp.npart;i++) P[i].sub = 0;

  for(i=0;i<nfrac;i++)
  {
    fprintf(stdout, "\nBegins Identification : Step %d of %d \n",i+1,nfrac);
    
    iden.r0  = fof[i];
    iden.r0 *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0;
    iden.step = i;
    iden.nobj = cp.npart;

    if(iden.r0 <= cp.soft)
    {
      fprintf(stdout,"cambia Linking length = %f \n",iden.r0);
      iden.r0 = cp.soft;
      i = nfrac;
    }

    fprintf(stdout,"Linking length = %f \n",iden.r0);

    identification();

    Write_Groups(fof[i]);

    free(Temp.head);
    free(Temp.npgrup);
    free(Temp.ll);
  }

  /************* TERMINO LA IDENTIFICACION ***************/

  free(P);
  grid_free();

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}

void Write_Groups(double fof)
{
  int i,j,k,id,npar,gn,save_sub;
  type_real dx,dy,dz,xc,yc,zc;
  char filename[200];
  FILE *pfout, *pfcentros;
  #ifdef FILE_ASCII
    FILE *pfcentros_ascii;
  #endif

  i = iden.ngrupos-1; // LE RESTO UNO POR EL GRUPO 0

  ///////////////////////////////////////////////////////
  sprintf(filename,"%.2d_%.2f_fof.bin",snap.num,fof);
  pfout=fopen(filename,"w");
  fwrite(&i,sizeof(int),1,pfout);
  //////////////////////////////////////////////////////
  sprintf(filename,"%.2d_%.2f_centros.bin",snap.num,fof);
  pfcentros=fopen(filename,"w");
  fwrite(&i,sizeof(int),1,pfcentros);
  //////////////////////////////////////////////////////
  #ifdef FILE_ASCII
    sprintf(filename,"%.2d_%.2f_centros.dat",snap.num,fof);
    pfcentros_ascii=fopen(filename,"w");
  #endif  
  //////////////////////////////////////////////////////

  npar = gn = 0;

  for(i=1;i<iden.ngrupos;i++)
  {

    j = 0;
    id = k = Temp.head[i];
    xc = yc = zc = 0.0;    

    if(iden.step==0)
      save_sub = i;
    else
      save_sub = P[k].sub;

    fwrite(&save_sub,sizeof(int),1,pfout);
    fwrite(&i,sizeof(int),1,pfout);
    fwrite(&Temp.npgrup[i],sizeof(int),1,pfout);

    while(k != -1)
    {

      // cuidado con el orden {pos[i]-centro} en este caso
      dx = P[k].Pos[0] - P[id].Pos[0];
      dy = P[k].Pos[1] - P[id].Pos[1];
      dz = P[k].Pos[2] - P[id].Pos[2];

      #ifdef PERIODIC
      dx = dx > cp.lbox*0.5 ? dx-cp.lbox : dx;
      dy = dy > cp.lbox*0.5 ? dy-cp.lbox : dy;
      dz = dz > cp.lbox*0.5 ? dz-cp.lbox : dz;
  
      dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;
      dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;
      dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
      #endif

      xc += dx;
      yc += dy;
      zc += dz;      

      P[k].sub = P[k].gr;
      fwrite(&P[k].id,sizeof(int),1,pfout);
      k = Temp.ll[k];
      j++;
    }
    
    #ifdef DEBUG
    assert(j == Temp.npgrup[i]);
    #endif

    xc /= (type_real)Temp.npgrup[i];
    yc /= (type_real)Temp.npgrup[i];
    zc /= (type_real)Temp.npgrup[i];

    xc += P[id].Pos[0];
    yc += P[id].Pos[1];
    zc += P[id].Pos[2];

    xc += pmin[0];
    yc += pmin[1];
    zc += pmin[2];
         
    #ifdef PERIODIC
    xc = xc<0 ? cp.lbox+(float)fmod(xc,cp.lbox) : (float)fmod(xc,cp.lbox);
    yc = yc<0 ? cp.lbox+(float)fmod(yc,cp.lbox) : (float)fmod(yc,cp.lbox);
    zc = zc<0 ? cp.lbox+(float)fmod(zc,cp.lbox) : (float)fmod(zc,cp.lbox);
    #endif

    fwrite(&save_sub,sizeof(int),1,pfcentros);
    fwrite(&i,sizeof(int),1,pfcentros);
    fwrite(&xc,sizeof(float),1,pfcentros);
    fwrite(&yc,sizeof(float),1,pfcentros);
    fwrite(&zc,sizeof(float),1,pfcentros);
    fwrite(&Temp.npgrup[i],sizeof(int),1,pfcentros);

    #ifdef FILE_ASCII
      fprintf(pfcentros_ascii,"%d %d %f %f %f %d\n",save_sub,i,xc,yc,zc,Temp.npgrup[i]);
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

  fprintf(stdout,"num de grupos %d num de particulas en grupos %d\n",gn,npar);
  fflush(stdout);

  return;
}
