#ifndef LIST_H
#define LIST_H

struct data_list
{
  type_real  r;
  type_int  id;
};

struct list
{
  struct data_list data;
  struct list *next;
};

extern void push_list(struct list **root, struct data_list dat);
extern void free_list(struct list *root);
extern void remove_duplicates(struct list *root);
extern void print_lista(struct list *lista);
extern void allocate_ids(struct list *root, type_int *array);
extern void count_list(struct list *root, type_int *c);

#endif
