#ifndef PROP_H
#define PROP_H

#ifndef VARIABLES_H
  #include "variables.h"
#endif

struct propiedades_st
{
  type_int   npart;
  type_real  *pos;
  type_real  pcm[3];
  type_real  aa, bb, cc;
  type_real  evec[3][3];
#ifdef STORE_VELOCITIES  
  type_real  *vel;
  type_real  vcm[3], sig[3];
  type_real  r200, m200, v200;
  type_real  rvir, mvir, vvir;
  type_real  vmax;
  type_real  L[3]; 
  type_real  aa_vel, bb_vel, cc_vel;
  type_real  evec_vel[3][3];
#ifdef COMPUTE_EP
  type_real  mostbound[3];
  type_real  Ep, Ec, lambda;
#endif
#endif
};

extern void propiedades(struct propiedades_st *Prop);

#endif
