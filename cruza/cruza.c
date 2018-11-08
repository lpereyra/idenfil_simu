#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "colores.h"

int main(int argc, char **argv)
{
  type_int c, i, j, k, *test;
  double start,end;
  char filename[200];
  FILE *pf;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  read_segment(0);
  read_grup_fof(0);

  read_segment(1);
  read_grup_fof(1);

  c = 0;
  test = (type_int *) malloc(cp.nseg_cut*sizeof(type_int));

  for(j=0;j<cp.nseg_cut;j++)
  {

    test[j] = cp.nseg;

    if(Seg_cut[j].flag != 2) continue;

    #pragma omp parallel for num_threads(NTHREADS) \
    schedule(dynamic) default(none) private(i) \
    shared(cp,j,Seg,Gr,Seg_cut,Gr_cut,test,stdout) reduction(+:c)
    for(i=0;i<cp.nseg;i++)
    {

      if(Seg[i].flag != 2) continue;

      if((Gr[Seg[i].list[0]].id == Gr_cut[Seg_cut[j].list[0]].id) && \
      (Gr[Seg[i].list[Seg[i].size-1]].id == Gr_cut[Seg_cut[j].list[Seg[j].size-1]].id))
      {

        //if(test[j]!=cp.nseg)
        //  fprintf(stdout,"%d %d %d -- %d %d %f %f %f %d -- %d %d %f %f %f %d\n", i, j, test[j], \
        //  Gr[Seg[i].list[0]].save,Gr[Seg[i].list[0]].id, \
        //  Gr[Seg[i].list[0]].Pos[0],Gr[Seg[i].list[0]].Pos[1],Gr[Seg[i].list[0]].Pos[2], \
        //  Gr[Seg[i].list[0]].NumPart, \
        //  Gr_cut[Seg_cut[j].list[0]].save,Gr_cut[Seg_cut[j].list[0]].id,\
        //  Gr_cut[Seg_cut[j].list[0]].Pos[0],Gr_cut[Seg_cut[j].list[0]].Pos[1],Gr_cut[Seg_cut[j].list[0]].Pos[2],\
        //  Gr_cut[Seg_cut[j].list[0]].NumPart);

        test[j] = i;
        c++;
      }
    }
  }
  
  fprintf(stdout,"coincidencias %d\n",c);
 
  sprintf(filename,"%s_coicidencias_%s.bin",fil.fname,fil_cut.fname);
  pf = fopen(filename,"w");

  k = c;
  fwrite(&c,sizeof(type_int),1,pf);        

  for(j=0;j<cp.nseg_cut;j++)
  {
    if(test[j] == cp.nseg) continue;
    fwrite(&test[j],sizeof(type_int),1,pf); // Num Filamento 
    fwrite(&j,sizeof(type_int),1,pf);       // Num Filamento Cut
    c--;
  }
  
  fprintf(stdout,"DEBE SER CERO %d\n",c); 
  
  if(c!=0)
  {
    c = k-c;
    rewind(pf);
    fwrite(&c,sizeof(type_int),1,pf);        
  }

  fclose(pf);
  
  //assert(c==0);

  for(i=0;i<cp.nseg;i++)
    free(Seg[i].list);
  free(Seg);
  free(Gr);

  for(i=0;i<cp.nseg_cut;i++)
    free(Seg_cut[i].list);
  free(Seg_cut);
  free(Gr_cut);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);

}
