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
type_real *d, type_real e, type_real f, type_real g, type_real h) 
{
  int i;
  struct node_sph *t = root;
  
  while(t->next != NULL)
      t = t->next;
  
  t->next = (struct node_sph *) malloc(sizeof(struct node_sph));
  t->next->a = a;
  t->next->b = b;
  for(i=0;i<3;i++)
  {
    t->next->c[i] = c[i];
    t->next->d[i] = d[i];
  } 
  t->next->e = e;
  t->next->f = f;
  t->next->g = g;
  t->next->h = h;


  t->next->next = NULL;
}

void push_array(struct node_sph *root, int a, type_real b, type_real *c, \
type_real *d, type_real e, type_real f, type_real g, type_real h) 
{
  int i;
  
  root->a = a;
  root->b = b;
  for(i=0;i<3;i++)
  {
    root->c[i] = c[i];
    root->d[i] = d[i];
  } 
  root->e = e;
  root->f = f;
  root->g = g;
  root->h = h;
  root->next = NULL;
}

/////////////////////////////////////////////////////////////////////////

#endif
