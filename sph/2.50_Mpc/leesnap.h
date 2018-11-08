#ifndef LEESNAP_H
#define LEESNAP_H

//Input and output files//
struct SnapST{
  int nfiles;
  char root[200], name[200]; 
  char box[200], seg[200]; 
  int num;
};

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
};

extern struct SnapST    snap;  

extern void read_segment(type_int NNN, type_real *fof);

#endif
