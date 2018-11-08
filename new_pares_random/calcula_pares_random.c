#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "variables.h"
#include "cosmoparam.h"
#include "colores.h"
#include "leesnap.h"
#include "calcula_pares_random.h"
#include <iterator>
#include <vector>
#include <set>

void crea_random(int NpartCut)
{
  type_int i,j,idnod,idnew;
  std::vector< std::vector<type_int> > idgrup;
  std::set<type_int>::iterator it;
  std::set<type_int> s;
  type_int *idvec;

  idvec = (type_int *) malloc(cp.ngrup*sizeof(type_int));

  for(i=0;i<cp.ngrup;i++)
    s.insert(Gr[i].NumPart);

  for(it = s.begin(); it != s.end(); it++)
    idgrup.push_back( std::vector<type_int> () );    

  for(i=0;i<cp.ngrup;i++)
  {
    it = s.find(Gr[i].NumPart);

    j = std::distance(s.begin(), it);

    idvec[i] = j;
    idgrup[j].push_back(i);
  }

  i = 0;
  for(j=0;j<s.size();j++)
    i += idgrup[j].size();

  s.empty(); // libero

  assert(i==cp.ngrup);

  srand(80); // guarda la semilla

  for(i=0;i<cp.nseg;i++)
  {

    idnod = Seg[i].list[0];
    j = idgrup[idvec[idnod]].size(); // size
    j = rand() % j;                  // tiro el random
    idnew = idgrup[idvec[idnod]][j]; // elijo otro halo random;
    assert(Gr[idnod].NumPart == Gr[idnew].NumPart);
    Seg[i].list[0] = idnew; // sobrescribo

    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    
    idnod = Seg[i].list[1];
    j = idgrup[idvec[idnod]].size(); // size
    j = rand() % j;                  // tiro el random
    idnew = idgrup[idvec[idnod]][j]; // elijo otro halo random;
    assert(Gr[idnod].NumPart == Gr[idnew].NumPart);
    Seg[i].list[1] = idnew; // sobrescribo

  }
 
  while(!idgrup.empty())
    idgrup.pop_back();
  free(idvec);

  return;
}


