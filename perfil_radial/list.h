#ifndef LIST_H
#define LIST_H

/////////////////////////////////////////////////////////////////////////

struct node_rad{
  int       a;
  type_real b;
  type_real c[3];
  type_real d[3];
  type_real e;
  type_real f;
  struct node_rad *next;
};

/////////////////////////////////////////////////////////////////////////

void push_list(struct node_rad *root, int a, type_real b, type_real *c, \
type_real *d, type_real e, type_real f) 
{
  int i;
  struct node_rad *t = root;
  
  while(t->next != NULL)
      t = t->next;
  
  t->next = (struct node_rad *) malloc(sizeof(struct node_rad));
  t->next->a = a;
  t->next->b = b;
  for(i=0;i<3;i++)
  {
    t->next->c[i] = c[i];
    t->next->d[i] = d[i];
  } 
  t->next->e = e;
  t->next->f = f;


  t->next->next = NULL;
}

void push_array(struct node_rad *root, int a, type_real b, type_real *c, \
type_real *d, type_real e, type_real f) 
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
  root->next = NULL;
}

/////////////////////////////////////////////////////////////////////////

#endif
