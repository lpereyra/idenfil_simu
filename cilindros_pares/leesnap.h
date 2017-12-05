#ifndef LEESNAP_H
#define LEESNAP_H

void leeheader(char *nombrefile);
void lee(char *filename, struct particle_data *Q, int *ind);
void read_gadget();
void read_pares(int NNN, type_real *fof);
void read_grup_fof(type_real *fof);

type_real pmin[3], pmax[3];

/*Input and output files*/
struct SnapST{
  int nfiles;
  char root[200], name[200];
  int num;
} snap;

struct io_header{
  int      npart[N_part_types];
  double   mass[N_part_types];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[N_part_types];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- N_part_types*4- N_part_types*8- 2*8- 2*4- N_part_types*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header;

#endif
