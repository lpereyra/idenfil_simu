#ifndef CALCULA_H
#define CALCULA_H

void crea_random(int NpartCut);
void sort_halos();
void reorder_grups(int *Id);
int compare_descend(const void *a, const void *b);
void find_cut(int *k, int NCut);

struct sort_prop{
  int index;
  type_real mass;
};

#endif
