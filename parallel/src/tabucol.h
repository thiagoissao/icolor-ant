/**********************************************************
 * Created: Seg 08 Ago 2011 17:22:32 BRT
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 **********************************************************
 *
 * Local search Tabucol [1987, Hertz]
 * Works with infeasible k-colorings.
 * Can be used with either a react or a dynamic tabu tenure.
 *
 ***********************************************************/

#ifndef __TABUCOL_H
#define __TABUCOL_H

#define TABUCOL_REACTIVE 0
#define TABUCOL_DYNAMIC 1

#define TABUCOL_CYCLES 1000000
#define TABUCOL_CONVERGENCE_CYCLES 10000
#define TABUCOL_CHANGE_SCHEME_ITERATIONS 100
#define TABUCOL_DIFF_SCHEME_ITERATIONS 1

struct tabucol_t {
  int tl_style;
  int cycles;
  int convergence_cycles;
  int change_scheme_iterations;
  int diff_scheme_iterations;
  double spent_time;
};

typedef struct tabucol_conflicts_t tabucol_conflicts_t;
struct tabucol_conflicts_t {
  /* Number of conflicts for each color and node: (k+1)x(n+1) */
  int **conflicts;
  /* Tabu status for each color and node: (k+1)x(n) */
  int **tabu_status;
  /* Nodes that have conflicts. Position 0 keeps quantity: (k+1)x(n+1) */
  int *nodes_in_conflict;
  /* Index of nodes in array nodes_in_conflicts */
  int *conf_position;
};

typedef struct tabucol_t tabucol_t;

tabucol_t *tabucol_info;

void tabucol_printbanner(void);
void tabucol_malloc(void);
void tabucol_show_solution(void);
void tabucol(tabucol_conflicts_t **tabucol_conflicts, gcp_solution_t *solution,
             int max_cycles, int type_of_tl);

#endif /* __TABUCOL_H */
