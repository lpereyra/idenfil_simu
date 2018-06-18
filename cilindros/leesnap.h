#ifndef LEESNAP_H
#define LEESNAP_H

extern void read_gadget(void);
extern void read_segment(const type_int NNN, const type_real *fof);
extern void read_grup_fof(const type_real *fof);
#ifdef VEL_RELATIVA
  extern void grupos_fof(const type_real *fof);
#endif

type_real pmin[3], pmax[3];

/*Input and output files*/
struct SnapST{
  int nfiles;
  char root[200], name[200];
  type_int num;
} snap;

struct io_header{
  type_int npart[N_part_types];
  double   mass[N_part_types];
  double   time;
  double   redshift;
  type_int flag_sfr;
  type_int flag_feedback;
  type_int npartTotal[N_part_types];
  type_int flag_cooling;
  type_int num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- N_part_types*4- N_part_types*8- 2*8- 2*4- N_part_types*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header;

#endif
