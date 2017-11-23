#ifndef IDENTIFICADOR_H
#define IDENTIFICADOR_H

struct iden_st
{
	double r0;        /*Linking Length para el FoF*/
	int    nobj;
	int    step;
  int    ngrupos;
} iden;

struct temporary 
{
	int *head;
	int *ll;
	unsigned int *npgrup;
} Temp;

void identification(void);
int Raiz(int i);
void linkedlist(int *grup_thread);
#ifdef LOCK
void Unir(int u, int v, omp_lock_t *lock);
void busv_rec(int i, int *test, omp_lock_t *lock);
#else
void Unir(int u, int v);
void busv_rec(int i, int *test);
#endif


#endif
