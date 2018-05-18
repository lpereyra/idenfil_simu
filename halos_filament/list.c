#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "variables.h"
#include "list.h"

void push_list(struct list **root, struct data_list dat)
{
  struct list * new = (struct list * ) malloc(sizeof(struct list));

  new->data = dat;
  new->next = (* root);
  (* root)  = new;

  /*
   *
  struct list * new = (struct list * ) malloc(sizeof(struct list));

  new->id   = id;
  new->next = NULL;

  if (*root == NULL)
  {
     *root = new;
     return;
  } 

  struct list *t = *root;

  while(t->next != NULL)
      t = t->next;
  
  t->next = new;
  */

  return;
}

void free_list(struct list *root)
{
  struct list *t;

  while(root != NULL) 
  {
      t = root;
      root = root->next;           
      free(t); // delete saved pointer.
  }
  
  return;

}

void remove_duplicates(struct list *root)
{
  struct list *p1, *p2, *sp;

  p1 = root;

  // ESTO ES N CUADRADO
  while (p1 != NULL && p1->next != NULL)
  {
      p2 = p1;
 
      while (p2->next != NULL)
      {
          if (p1->data.id == p2->next->data.id)
          { 

            // swap para quedarse con la menor distancia
            // guarda aca!!!!
            if(p1->data.r > p2->next->data.r)
              p1->data = p2->next->data;

            sp = p2->next;
            p2->next = p2->next->next;
            free(sp);
          }else{
            p2 = p2->next;
          }
      }
      p1 = p1->next;
  }

  return;
}

void print_lista(struct list *lista)
{
  while(lista != NULL)
  {
      fprintf(stdout,"(%d,%f),",lista->data.id, lista->data.r);
      lista = lista->next;
  }
  fprintf(stdout,"\n");
}

void count_list(struct list *root, type_int *c)
{
  *c = 0;
  while(root != NULL) 
  {
      root = root->next;          
      *c = *c + 1;
  }
  
  return;
}

void alocate_list(struct list *root, struct data_list *array)
{
  assert(array != NULL);

  type_int c = 0;

  while(root != NULL) 
  {
    array[c] = root->data;
    c = c + 1;
    root = root->next;          
  }
  
  return;
}

