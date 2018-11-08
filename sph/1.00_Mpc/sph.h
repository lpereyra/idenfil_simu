#ifndef SPH_H
#define SPH_H

extern void change_positions();
extern void init_sph(const type_int i);

extern "C" {
 void arvo_(int*,double*,double*,double*,double*,double*);
}

#endif
