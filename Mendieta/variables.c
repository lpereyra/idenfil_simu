#include "variables.h"

struct cosmoparam cp;
struct SnapST snap;
struct gridst grid;
#ifdef COLUMN
  struct particle_data P;
#else
  struct particle_data *P;
#endif
type_int  nfrac;
type_real *fof;
type_real pmin[3], pmax[3];
char message[200];
