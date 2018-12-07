#ifndef CHULL_H
#define CHULL_H

#include "variables.h"

#define SWAP(t,x,y)	{ t = x; x = y; y = t; }

//char *malloc();

#define NEW(p,type)	if ((p=(type *) malloc (sizeof(type))) == NULL) {\
				printf ("Out of Memory!\n");\
				exit(0);\
			}

#define FREE(p)		if (p) { free ((char *) p); p = NULL; }


#define ADD( head, p )  if ( head )  { \
				p->next = head; \
				p->prev = head->prev; \
				head->prev = p; \
				p->prev->next = p; \
			} \
			else { \
				head = p; \
				head->next = head->prev = p; \
			}

#define DELETE( head, p ) if ( head )  { \
				if ( head == head->next ) \
					head = NULL;  \
				else if ( p == head ) \
					head = head->next; \
				p->next->prev = p->prev;  \
				p->prev->next = p->next;  \
				FREE( p ); \
			} 

/*Define Boolean type */
typedef	enum { FALSE, TRUE }	bool;

/* Define vertex indices. */
#define X 0
#define Y 1
#define Z 2

/* Define structures for vertices, edges and faces */
typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;

typedef struct tEdgeStructure tsEdge;
typedef tsEdge *tEdge;

typedef struct tFaceStructure tsFace;
typedef tsFace *tFace;

struct tVertexStructure {
   type_int  idx;
   type_int	 vnum;
   tEdge     duplicate;	        /* pointer to incident cone edge (or NULL) */
   bool      onhull;		        /* T iff point on hull. */
   bool	     mark;	            /* T iff point already processed. */
   tVertex   next, prev;
};

struct tEdgeStructure {
   tFace    adjface[2];
   tVertex  endpts[2];
   tFace    newface;            /* pointer to incident cone face. */
   bool     delete;		/* T iff edge should be delete. */
   tEdge    next, prev;
};

struct tFaceStructure {
   tEdge    edge[3];
   tVertex  vertex[3];
   bool	    visible;	        /* T iff face visible from new point. */
   tFace    next, prev;
};

/* Define flags */
#define ONHULL   	TRUE
#define REMOVED  	TRUE
#define VISIBLE  	TRUE
#define PROCESSED	TRUE

extern void build_chull();

#endif
