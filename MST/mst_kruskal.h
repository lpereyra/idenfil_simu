#ifndef KRUSKAL 
#define KRUSKAL 

int  Root(int i, int *Padre);
void DLU(int id, int pre, std::vector<std::vector<int> > &adj, int *Padre, int *Rank);
void DL(std::vector<std::vector<int> > &vec, int *Padre, int *Rank);
bool DFS(int i, std::vector<std::vector<int> > &vec, int cut);
void Delete_Branch(int i, std::vector<std::vector<int> > &vec);
void Podado(int level, std::vector<std::vector<int> > &vec);
void Kruskal(int *Padre, int *Rank, std::vector<std::pair<type_real,std::pair<int,int> > > &edges, std::vector<std::vector<int> > &adjacency_list);

#endif
