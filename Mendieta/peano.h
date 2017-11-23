#ifndef PEANO_HILBERT 
#define PEANO_HILBERT 

typedef  long long  peanokey;    
#define  BITS_PER_DIMENSION 18	

struct peano_hilbert_data
{
  peanokey key;
  int index;
};

void      peano_hilbert(void);
peanokey  peano_hilbert_key(int x, int y, int z, int bits);
int       compare_key(const void *a, const void *b);
void reorder_particles(int *Id);

#endif
