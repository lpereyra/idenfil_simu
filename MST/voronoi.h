#ifndef VORONOI 
#define VORONOI 

bool dist_segment(double fof, int idg, type_real x, type_real y, type_real z);
void Voronoi_Grupos(const type_real fof, std::vector<std::pair<type_real,std::pair<int,int> > > &edges);
void init_edges(std::vector<std::pair<int,int> > &pares);
int inside_fof(const type_real fof, const type_int id, const type_real rcil, const type_real rcil2, const type_real r, const type_real vdir[]);
int control(const type_real fof, const type_real rlong_2, const type_real rcil2, const type_real rsph, const int idg, const type_real versor[], const type_real Pos_cent[]);

#endif

