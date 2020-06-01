#ifndef OCTREE_H
#define OCTREE_H

#ifndef VARIABLES_H
  #include "variables.h"
#endif

#ifndef PROP_H
  #include "propiedades.h"
#endif

struct NODE{ 
  type_real s[3];                     /* center of mass */
  long  partind;
  type_real center[3],len;            /* center and sidelength of treecubes */
  float mass;                         /* mass*/
  float oc;                           /* variable for opening criter*/
  struct NODE *next,*sibling,*father,*suns[8];
};

extern void force_treeallocate(int maxnodes);
extern void force_treefree(void);

extern int force_treebuild(struct propiedades_st *Prop, float thetamax);
extern type_real force_treeevaluate_potential(const type_real pos[]);

#endif
