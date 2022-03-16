/**********************************************************
 * Created: Seg 10 Out 2011 17:23:49 BRT
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 * ANT_FIXED_K
 * * Constructive method for each ant in an ACO algorithm for k-GCP
 *
 **********************************************************/

#ifndef __ANT_FIXED_K_H
#define __ANT_FIXED_K_H

typedef struct ant_fixed_k_t ant_fixed_k_t;
struct ant_fixed_k_t {
  /* ANT_FIXED_K data */
  double *probb;
  double **trail;
  float alpha;
  float beta;
};

struct ant_fixed_k_t *afk_initialize_data(float p_alpha, float p_beta);
void ant_fixed_k(ant_fixed_k_t **ant_fixed_k, gcp_solution_t *solution,
                 double **pheromone);

#endif /* __ANT_FIXED_K_H */
