#ifndef _FAST3TREE_C_
#define _FAST3TREE_C_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>

#ifndef FAST3TREE_DIM
  #define FAST3TREE_DIM 3
#endif /* FAST3TREE_DIM */

#ifndef POINTS_PER_LEAF
  #define POINTS_PER_LEAF 40
#endif  /* POINTS_PER_LEAF */

#ifdef FAST3TREE_FLOATTYPE
  #define float FAST3TREE_FLOATTYPE
#endif /* FAST3TREE_FLOATTYPE */

#ifdef PERIODIC
  #define FAST3TREE_MARKED 1 //Flag for marked nodes
#endif

struct tree3_node {
  float min[FAST3TREE_DIM], max[FAST3TREE_DIM];
  int64_t num_points;
  int16_t div_dim, flags;
  struct tree3_node *left, *right, *parent;
  struct grup_data *points;
};

struct fast3tree {
  struct grup_data *points;
  int64_t num_points;
  struct tree3_node *root;
  int64_t num_nodes;
  int64_t allocated_nodes;
};

struct fast3tree_results {
  int64_t num_points;
  int64_t num_allocated_points;
  struct grup_data **points;
};

#undef fast3tree_init
#define fast3tree_init init
struct fast3tree *fast3tree_init(int64_t n, struct grup_data *p);

#undef fast3tree_rebuild
#define fast3tree_rebuild rebuild
void fast3tree_rebuild(struct fast3tree *t, int64_t n, struct grup_data *p);

#undef fast3tree_results_init
#define fast3tree_results_init results_init
struct fast3tree_results *fast3tree_results_init(void);

#undef fast3tree_find_sphere
#define fast3tree_find_sphere find_sphere
static inline void fast3tree_find_sphere(struct fast3tree *t,
			   struct fast3tree_results *res, float c[FAST3TREE_DIM], float r);

#ifdef PERIODIC

  #undef fast3tree_find_sphere_periodic
  #define fast3tree_find_sphere_periodic find_sphere_periodic
  int fast3tree_find_sphere_periodic(struct fast3tree *t,
			   struct fast3tree_results *res, float c[FAST3TREE_DIM], float r);

#endif

#undef fast3tree_results_clear
#define fast3tree_results_clear results_clear
void fast3tree_results_clear(struct fast3tree_results *res);

#undef fast3tree_results_free
#define fast3tree_results_free results_free
void fast3tree_results_free(struct fast3tree_results *res);

#ifndef __APPLE__
  #ifndef isfinite
    #define isfinite(x) finitef(x)
  #endif
#endif

/* PRIVATE METHODS */
#undef _fast3tree_check_realloc
#define _fast3tree_check_realloc _fast3tree_check_realloc
void *_fast3tree_check_realloc(void *ptr, size_t size, char *reason);

#undef _fast3tree_build
#define _fast3tree_build _fast3tree_build
void _fast3tree_build(struct fast3tree *t);

struct fast3tree *fast3tree_init(int64_t n, struct grup_data *p) {
  struct fast3tree *new=NULL;
  new = _fast3tree_check_realloc(new,sizeof(struct fast3tree), "Allocating fast3tree.");
  if (!new) return 0;
  memset(new, 0, sizeof(struct fast3tree));
  fast3tree_rebuild(new, n, p);
  return new;
}

void fast3tree_rebuild(struct fast3tree *t, int64_t n, struct grup_data *p) {
  t->points = p;
  t->num_points = n;
  _fast3tree_build(t);
}

#undef fast3tree_free
#define fast3tree_free fast3tree_free
void fast3tree_free(struct fast3tree **t) {
  if (!t) return;
  struct fast3tree *u = *t;
  if (u) {
    free(u->root);
    free(u);
  }
  *t = NULL;
}

/* Accurate, fast */
#undef _fast3tree_box_not_intersect_sphere
#define _fast3tree_box_not_intersect_sphere _fast3tree_box_not_intersect_sphere
static inline int _fast3tree_box_not_intersect_sphere(const struct tree3_node *node, const float c[FAST3TREE_DIM], const float r) {
  int i;
  float d = 0, e;
  const float r2 = r*r;
  for (i=0; i<FAST3TREE_DIM; i++) {
    if ((e = c[i]-node->min[i])<0) {
      d+=e*e;
      if (d >= r2) return 1;
    } else if ((e = c[i]-node->max[i])>0) {
      d+=e*e;
      if (d >= r2) return 1;
    }
  }
  return 0;
}

/* Fast, accurate. */
#undef _fast3tree_box_inside_sphere
#define _fast3tree_box_inside_sphere _fast3tree_box_inside_sphere
inline int _fast3tree_box_inside_sphere(const struct tree3_node *node, const float c[FAST3TREE_DIM], const float r) {
  int i;
  float dx, dx2, dist = 0, r2 = r*r;
  if (fabs(c[0]-node->min[0]) > r) return 0; //Rapid short-circuit.
  for (i=0; i<FAST3TREE_DIM; i++) {
    dx = node->min[i] - c[i];
    dx *= dx;
    dx2 = c[i]-node->max[i];
    dx2 *= dx2;
    if (dx2 > dx) dx = dx2;
    dist += dx;
    if (dist > r2) return 0;
  }
  return 1;
}

#ifdef PERIODIC

  #undef _fast3tree_sphere_inside_box
  #define _fast3tree_sphere_inside_box _fast3tree_sphere_inside_box
  static inline int _fast3tree_sphere_inside_box(const struct tree3_node *node, const float c[FAST3TREE_DIM], const float r) {
    int i;
    for (i=0; i<FAST3TREE_DIM; i++) {
      if (c[i]-r < node->min[i]) return 0;
      if (c[i]+r > node->max[i]) return 0;
    }
    return 1;
  }

#endif

#undef _fast3tree_check_results_space
#define _fast3tree_check_results_space _fast3tree_check_results_space
static inline void _fast3tree_check_results_space(const struct tree3_node *n, struct fast3tree_results *res) {
  if (res->num_points + n->num_points > res->num_allocated_points) {
    res->num_allocated_points = res->num_points + n->num_points + 1000;

    res->points = _fast3tree_check_realloc(res->points, 
    res->num_allocated_points * sizeof(struct grup_data *), "Allocating fast3tree results");
  }
}

#undef _fast3tree_find_sphere 
#define _fast3tree_find_sphere _fast3tree_find_sphere
void _fast3tree_find_sphere(const struct tree3_node *n, struct fast3tree_results *res, float c[FAST3TREE_DIM], const float r) {
  int64_t i,j;
  float r2, dist, dx;

  if (_fast3tree_box_not_intersect_sphere(n,c,r)) return;
#if FAST3TREE_DIM < 6
  if (_fast3tree_box_inside_sphere(n,c,r)) { /* Entirely inside sphere */
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++)
      res->points[res->num_points+i] = n->points+i;
    res->num_points += n->num_points;
    return;
  }
#endif /* FAST3TREE_DIM < 6 */

  if (n->div_dim < 0) { /* Leaf node */
    r2 = r*r;
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++) {
      j = dist = 0;
      float *Pos = n->points[i].Pos;
      for (; j<FAST3TREE_DIM; j++) {
	dx = c[j]-Pos[j];
	dist += dx*dx;
      }
      if (dist < r2) {
	res->points[res->num_points] = n->points + i;
	res->num_points++;
      }
    }
    return;
  }
  _fast3tree_find_sphere(n->left, res, c, r);
  _fast3tree_find_sphere(n->right, res, c, r);
}

struct fast3tree_results *fast3tree_results_init(void) {
  struct fast3tree_results *res = 
    _fast3tree_check_realloc(NULL, sizeof(struct fast3tree_results), "Allocating fast3tree results structure.");
  res->points = NULL;
  res->num_points = 0;
  res->num_allocated_points = 0;
  return res;
}

static inline void fast3tree_find_sphere(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r) {
  res->num_points = 0;
  _fast3tree_find_sphere(t->root, res, c, r);
}

#ifdef PERIODIC

  /* Guaranteed to be stable with respect to floating point round-off errors.*/
  #undef _fast3tree_find_sphere_offset
  #define _fast3tree_find_sphere_offset _fast3tree_find_sphere_offset
  void _fast3tree_find_sphere_offset(struct tree3_node *n, struct fast3tree_results *res, \
                                     float c[FAST3TREE_DIM], float c2[FAST3TREE_DIM], float o[FAST3TREE_DIM], \
                                     const float r, const int marked, const int do_marking) {
    int64_t i,j;
    float r2, dist, dx;
  
    int64_t onlyone = (marked && (n->flags & FAST3TREE_MARKED)) ? 1 : 0;
  
    if (_fast3tree_box_not_intersect_sphere(n,c2,r*1.01)) return;
    if (_fast3tree_box_inside_sphere(n,c2,r*0.99)) { /* Entirely inside sphere */
      if (do_marking) n->flags |= FAST3TREE_MARKED;
      _fast3tree_check_results_space(n,res);
      if (onlyone) {
        res->points[res->num_points++] = n->points;
        return;
      }
      for (i=0; i<n->num_points; i++)
        res->points[res->num_points+i] = n->points + i;
      res->num_points += n->num_points;
      return;
    }
  
    if (n->div_dim < 0) { /* Leaf node */
      r2 = r*r;
      _fast3tree_check_results_space(n,res);
      for (i=0; i<n->num_points; i++) {
        j = dist = 0;
        float *Pos = n->points[i].Pos;
        for (; j<FAST3TREE_DIM; j++) {
  	      dx = o[j] - fabs(c[j]-Pos[j]);
  	      dist += dx*dx;
        }
        if (dist < r2) {
  	      res->points[res->num_points] = n->points + i;
  	      res->num_points++;
  	      if (onlyone) return;
        }
      }
      return;
    }
    int64_t cur_points = res->num_points;
    _fast3tree_find_sphere_offset(n->left, res, c, c2, o, r, marked, do_marking);
    if (onlyone && (cur_points < res->num_points)) return;
    _fast3tree_find_sphere_offset(n->right, res, c, c2, o, r, marked, do_marking);
  }

  #undef _fast3tree_find_sphere_periodic_dim
  #define _fast3tree_find_sphere_periodic_dim _fast3tree_find_sphere_periodic_dim
  void _fast3tree_find_sphere_periodic_dim(struct fast3tree *t, struct fast3tree_results *res, \
                                           float c[FAST3TREE_DIM], float c2[FAST3TREE_DIM], float o[FAST3TREE_DIM], \
                                           float r, float dims[FAST3TREE_DIM], int dim, int marked, int do_marking) {
    float c3[FAST3TREE_DIM];
    if (dim<0) {
      _fast3tree_find_sphere_offset(t->root, res, c, c2, o, r, marked, do_marking);
      return;
    }
    memcpy(c3, c2, sizeof(float)*FAST3TREE_DIM);
    o[dim]=0;
    _fast3tree_find_sphere_periodic_dim(t, res, c, c3, o, r, dims, dim-1, marked, do_marking);
    if (c[dim]+r > t->root->max[dim]) {
      c3[dim] = c[dim]-dims[dim];
      o[dim] = dims[dim];
      _fast3tree_find_sphere_periodic_dim(t, res, c, c3, o, r, dims, dim-1, marked, do_marking);
    }
    if (c[dim]-r < t->root->min[dim]) {
      c3[dim] = c[dim]+dims[dim];
      o[dim] = dims[dim];
      _fast3tree_find_sphere_periodic_dim(t, res, c, c3, o, r, dims, dim-1, marked, do_marking);
    }
  }

  int fast3tree_find_sphere_periodic(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r) {
    float dims[FAST3TREE_DIM], o[FAST3TREE_DIM];
    int i;
    
    if (_fast3tree_sphere_inside_box(t->root, c, r)) {
      fast3tree_find_sphere(t, res, c, r);
      return 2;
    }
  
    for (i=0; i<FAST3TREE_DIM; i++) {
      dims[i] = t->root->max[i] - t->root->min[i];
      o[i] = 0;
      if (r*2.0 > dims[i]) return 0; //Avoid wraparound intersections.
    }
  
    res->num_points = 0;
    _fast3tree_find_sphere_periodic_dim(t, res, c, c, o, r, dims, FAST3TREE_DIM-1, 0, 0);
    return 1;
  }

#endif

void fast3tree_results_clear(struct fast3tree_results *res) {
  if (res->points) free(res->points);
  memset(res, 0, sizeof(struct fast3tree_results));
}

void fast3tree_results_free(struct fast3tree_results *res) {
  if (!res) return;
  if (res->points) free(res->points);
  memset(res, 0, sizeof(struct fast3tree_results));
  free(res);
}

#undef _fast3tree_find_largest_dim
#define _fast3tree_find_largest_dim _fast3tree_find_largest_dim
static inline int64_t _fast3tree_find_largest_dim(float * min, float * max) {
  int64_t i, dim = FAST3TREE_DIM-1;
  float d = max[FAST3TREE_DIM-1]-min[FAST3TREE_DIM-1], d2;
  for (i=0; i<(FAST3TREE_DIM-1); i++) {
    d2 = max[i] - min[i];
    if (d2 > d) { d=d2; dim = i; }
  }
  return dim;
}

#undef _fast3tree_sort_dim_pos
#define _fast3tree_sort_dim_pos _fast3tree_sort_dim_pos
static inline int64_t _fast3tree_sort_dim_pos(struct tree3_node * node,
					      float balance) {
  int64_t dim = node->div_dim = 
    _fast3tree_find_largest_dim(node->min, node->max);
  struct grup_data *p = node->points;
  int64_t i,j=node->num_points-1;
  struct grup_data temp;
  float lim = 0.5*(node->max[dim]+node->min[dim]);

  if (node->max[dim]==node->min[dim]) return node->num_points;

  for (i=0; i<j; i++) {
    if (p[i].Pos[dim] > lim) {
      temp = p[j];
      p[j] = p[i];
      p[i] = temp;
      i--;
      j--;
    }
  }

  if ((i==j) && (p[i].Pos[dim] <= lim)) i++;
  return i;
}

#undef _fast3tree_find_minmax
#define _fast3tree_find_minmax _fast3tree_find_minmax
static inline void _fast3tree_find_minmax(struct tree3_node *node) {
  int64_t i, j;
  float x;
  struct grup_data * p = node->points;
  assert(node->num_points > 0);
  for (j=0; j<FAST3TREE_DIM; j++) node->min[j] = node->max[j] = p[0].Pos[j];
  for (i=1; i<node->num_points; i++)  {
    for (j=0; j<FAST3TREE_DIM; j++) {
      x = p[i].Pos[j];
      if (x<node->min[j]) node->min[j] = x;
      else if (x>node->max[j]) node->max[j] = x;
    }
  }
}

#undef _fast3tree_split_node
#define _fast3tree_split_node _fast3tree_split_node
void _fast3tree_split_node(struct fast3tree *t, struct tree3_node *node) {
  int64_t num_left;
  struct tree3_node *null_ptr = NULL;
  struct tree3_node *left, *right;
  int64_t left_index, node_index;

  num_left = _fast3tree_sort_dim_pos(node, 1.0);
  if (num_left == node->num_points || num_left == 0) 
  { //In case all node points are at same spot
    node->div_dim = -1;
    return;
  }

  node_index = node - t->root;
  if ((t->num_nodes+2) > t->allocated_nodes) {
    t->allocated_nodes = t->allocated_nodes*1.05 + 1000;
    t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->allocated_nodes), "Tree nodes");
    node = t->root + node_index;
  }

  left_index = t->num_nodes;
  t->num_nodes+=2;

  node->left = null_ptr + left_index;
  node->right = null_ptr + (left_index + 1);

  left = t->root + left_index;
  right = t->root + (left_index + 1);
  memset(left, 0, sizeof(struct tree3_node)*2);

  right->parent = left->parent = null_ptr + node_index;
  left->num_points = num_left;
  right->num_points = node->num_points - num_left;
  left->div_dim = right->div_dim = -1;
  left->points = node->points;
  right->points = node->points + num_left;

  _fast3tree_find_minmax(left);
  _fast3tree_find_minmax(right);

  if (left->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, left);

  right = t->root + (left_index + 1);
  if (right->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, right);
}

#undef _fast3tree_rebuild_pointers
#define _fast3tree_rebuild_pointers _fast3tree_rebuild_pointers
void _fast3tree_rebuild_pointers(struct fast3tree *t) {
  int64_t i;
  struct tree3_node *nullptr = NULL;
  for (i=0; i<t->num_nodes; i++) {
#define REBUILD(x) x = t->root + (x - nullptr)
    REBUILD(t->root[i].left);
    REBUILD(t->root[i].right);
    REBUILD(t->root[i].parent);
#undef REBUILD
  }
}

#undef _fast3tree_build
#define _fast3tree_build _fast3tree_build
void _fast3tree_build(struct fast3tree *t) {
  int64_t i, j;
  struct tree3_node *root;
  struct grup_data tmp;
  t->allocated_nodes = (3+t->num_points/(POINTS_PER_LEAF/2));
  t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->allocated_nodes), "Tree nodes"); //Estimate memory load
  t->num_nodes = 1;

  //Get rid of NaNs / infs
  for (i=0; i<t->num_points; i++) {
    for (j=0; j<FAST3TREE_DIM; j++) if (!isfinite(t->points[i].Pos[j])) break;
    if (j<FAST3TREE_DIM) {
      tmp = t->points[i];
      t->num_points--;
      t->points[i] = t->points[t->num_points];
      t->points[t->num_points] = tmp;
      i--;
    }
  }

  root = t->root;
  memset(root, 0, sizeof(struct tree3_node));
  root->num_points = t->num_points;
  root->points = t->points;
  root->div_dim = -1;
  if (t->num_points) _fast3tree_find_minmax(root);
  for (j=0; j<FAST3TREE_DIM; j++) assert(isfinite(root->min[j]));
  for (j=0; j<FAST3TREE_DIM; j++) assert(isfinite(root->max[j]));

  if (root->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, root);

  t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->num_nodes), "Tree nodes");
  t->allocated_nodes = t->num_nodes;
  _fast3tree_rebuild_pointers(t);
}

#undef _fast3tree_check_realloc
#define _fast3tree_check_realloc _fast3tree_check_realloc
void *_fast3tree_check_realloc(void *ptr, size_t size, char *reason) {
  void *res = realloc(ptr, size);
  if ((res == NULL) && (size > 0)) {
    fprintf(stderr, "[Error] Failed to allocate memory (%s)!\n", reason);
    exit(1);
  }
  return res;
}

#undef _fast3tree_set_minmax
#define _fast3tree_set_minmax fast3tree_set_minmax
void _fast3tree_set_minmax(struct fast3tree *t, float min, float max) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    t->root->min[i] = min;
    t->root->max[i] = max;
  }
}

#undef float

#endif
