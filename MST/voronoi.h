#ifndef VORONOI 
#define VORONOI 

bool dist_segment(double fof, int idg, type_real x, type_real y, type_real z);
void Voronoi_Grupos(type_real fof, std::vector<std::pair<type_real,std::pair<int,int> > > &edges);

#endif

