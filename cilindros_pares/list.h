#ifndef LIST_H
#define LIST_H

/////////////////////////////////////////////////////////////////////////

struct node_sph{
  int       a;
  type_real b;
  type_real c[3];
  type_real d[3];
  type_real e;
  type_real f;
  type_real g;
  type_real h;
  struct node_sph *next;
};

/////////////////////////////////////////////////////////////////////////

void push_list(struct node_sph *root, int a, type_real b, type_real *c, \
type_real *d, type_real e, type_real f, type_real g, type_real h); 

void push_array(struct node_sph *root, int a, type_real b, type_real *c, \
type_real *d, type_real e, type_real f, type_real g, type_real h);
/////////////////////////////////////////////////////////////////////////

#endif
