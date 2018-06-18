#ifndef VORONOI 
#define VORONOI 

void Voronoi_Grupos(const type_real fof, std::vector<std::pair<type_real,std::pair<type_int,type_int> > > &edges);
void init_edges(std::vector<std::pair<type_int,type_int> > &pares);
//type_int inside_fof(const type_real fof, const type_int id, const type_real rcil, const type_real rcil2, const type_real r, const type_real vdir[]);
type_int inside_fof(const type_int id, const type_real rcil, const type_real rcil2, const type_real r, const type_real vdir[]);
//type_int control(const type_real fof, const type_real rlong_2, const type_real rcil2, const type_real rsph, const type_int idg, const type_real versor[], const type_real Pos_cent[]);

type_int control(const type_real rsph, const type_real rsph2, const type_int idg, const type_real Pos_cent[]);
#endif

