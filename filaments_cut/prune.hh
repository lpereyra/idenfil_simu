#ifndef PRUNE 
#define PRUNE 

#ifdef PRUNED
  extern void Podado(const type_int level, std::vector<std::vector<type_int> > &vec);
#endif
extern void DL(std::vector<std::pair<type_int,type_int> > &mass, std::vector<std::vector<type_int> > &vec, \
type_int * __restrict__ Padre, type_int * __restrict__ Rank);

#endif
