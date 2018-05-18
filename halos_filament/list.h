#ifndef LIST_H
#define LIST_H

#include "variables.h"

struct list
{
  struct data_list data;
  struct list *next;
};

void push_list(struct list **root, struct data_list dat);
void free_list(struct list *root);
void remove_duplicates(struct list *root);
void print_lista(struct list *lista);
void alocate_list(struct list *root, struct data_list *array);
void count_list(struct list *root, type_int *c);

#endif
