#ifndef KRUSKAL 
#define KRUSKAL 

type_int Root(type_int i, type_int *Padre);
void DLU(const type_int id, const type_int pre, std::vector<std::vector<type_int> > &adj, type_int *Padre, type_int *Rank);
void DL(std::vector<std::vector<type_int> > &vec, type_int *Padre, type_int *Rank);
type_int DFS(type_int i, std::vector<std::vector<type_int> > &vec, type_int cut);
void Delete_Branch(const type_int i, std::vector<std::vector<type_int> > &vec);
void Podado(type_int level, std::vector<std::vector<type_int> > &vec);
void Kruskal(type_int *Padre, type_int *Rank, std::vector<std::pair<type_real,std::pair<type_int,type_int> > > &edges, std::vector<std::vector<type_int> > &adjacency_list);

#endif
