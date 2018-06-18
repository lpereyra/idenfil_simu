#ifndef LIST_H
#define LIST_H

/////////////////////////////////////////////////////////////////////////

struct node_sph{
  type_int       a;
  type_real b;
  type_real c[3];
  type_real d[3];
  type_real e;
  type_real f;
  type_real g;
  type_real h;
  #ifdef LINKED_LIST
  struct node_sph *next;
  #endif
};

struct dat_struct
{
	type_real * nodo;
  type_real * rbin_centre;
	type_int **npart;
  type_real **rho_delta;
	type_real **mean;
  type_real **rms;
	type_real **mean_par;
  type_real **rms_par;
	type_real **mean_perp;
  type_real **rms_perp;
} *dat;

/////////////////////////////////////////////////////////////////////////
#ifdef LINKED_LIST
extern  void push_list(struct node_sph *root, const type_int a, const type_real b, const type_real * restrict c, \
const type_real * restrict d, const type_real e, const type_real f, const type_real g, const type_real h);
#endif

extern void push_array(struct node_sph *root, const type_int a, const type_real b, const type_real * restrict c, \
const type_real * restrict d, const type_real e, const type_real f, const type_real g, const type_real h); 
/////////////////////////////////////////////////////////////////////////

#endif
