#include <stdlib.h>
#include <stdio.h>

#include "variables.h"
#include "list.h"

#ifdef LINKED_LIST

  void push_list(struct node_sph *root, const type_int a, const type_real b, const type_real * restrict c, \
  const type_real * restrict d, const type_real e, const type_real f, const type_real g, const type_real h) 
  {
    type_int i;
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

#endif

void push_array(struct node_sph *root, const type_int a, const type_real b, const type_real * restrict c, \
const type_real * restrict d, const type_real e, const type_real f, const type_real g, const type_real h) 
{
  type_int i;
  
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
}
