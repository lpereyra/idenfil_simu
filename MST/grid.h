#ifndef GRID_H
#define GRID_H

#ifndef NGRIDMAX
#define NGRIDMAX 512
#endif

struct gridst
{
  int step;
	unsigned long ngrid;
	unsigned long nobj;
	int *llirst;
	int *ll;
};

extern struct gridst grid;

void grid_init(void);
void grid_build(void);
void grid_free(void);

#endif
