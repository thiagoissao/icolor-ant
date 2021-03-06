/*************************************************************************
 * Created: Qua 27 Jan 2011 15:17:39 BRST
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 * Modified By: Thiago Issao Yasunaka, thiagoyasunaka@hotmail.com
 * Mon 14 Fev 2022 22:13:25 BRST
 *************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>

#include "icolorant.h"
#include "../tabucol.h"
#include "../ant_fixed_k.h"
#include "../util.h"
#include "../tabucol.h"
#include "icolorant.utils.h"

void colorant_malloc(void) {

  aco_info = malloc_(sizeof(aco_t));
  aco_info->alpha = COLORANT_ALPHA;
  aco_info->beta = COLORANT_BETA;
  aco_info->rho = COLORANT_RHO;
  aco_info->ants = COLORANT_ANTS;
  aco_info->ratio = 0;
  aco_info->pheromone_scheme = PHEROMONE_SCHEME_1;
  aco_info->change_phero_scheme_iterations =
      COLORANT_CHANGE_PHERO_SCHEME_ITERATIONS;
  aco_info->memory_size = COLORANT_MEMORY_SIZE;
  aco_info->delta = COLORANT_DELTA;
  aco_info->gap = COLORANT_GAP;
  aco_info->gamma = COLORANT_GAMMA;
  aco_info->omega = COLORANT_OMEGA;
  aco_info->x = COLORANT_X;
  aco_info->y = COLORANT_Y;
  aco_info->iterations_alpha_beta = COLORANT_ITERATIONS_ALPHA_BETA;
}

void colorant_initialization(void) {

  aco_info->alpha_base = aco_info->alpha;
  aco_info->beta_base = aco_info->beta;

  if (get_flag(problem->flags, FLAG_ANTS_RATIO)) {
    aco_info->ratio = aco_info->ants;
    aco_info->ants = (problem->nof_vertices * aco_info->ants) / 100;
  }

  if (get_flag(problem->flags, FLAG_MEMORY_RATIO)) {
    aco_info->memory_ratio = aco_info->memory_size;
    aco_info->memory_size = (aco_info->memory_size * aco_info->ants) / 100;
    aco_info->memory_size =
        aco_info->memory_size < 1 ? 1 : aco_info->memory_size;
  }
}

int memory_length(aco_memory_t *memory) {
  aco_memory_t *item = memory;
  int length = 0;
  for (; item != NULL; item = item->tail, length++)
    ;
  return length;
}

void insert_into_memory(gcp_solution_t *solution,
                        gcp_solution_t *local_ant_memory_insert,
                        gcp_solution_t *local_ant_memory_remove,
                        aco_memory_t *memory) {
  aco_memory_t *item = malloc_(sizeof(aco_memory_t *));
  aco_memory_t *last, *previous;

  item->head = solution;
  item->tail = memory;
  memory = item;

  cpy_solution(solution, local_ant_memory_insert);

  if (aco_info->memory_size < memory_length(memory)) {
    last = memory->tail;
    previous = memory;
    for (; last && last->tail != NULL;
         last = last->tail, previous = previous->tail)
      ;
    previous->tail = NULL;
    cpy_solution(last->head, local_ant_memory_remove);
    set_flag(problem->flags, FLAG_MEMORY_REMOVE);
    free(last);
  }

  print_memory(memory, local_ant_memory_remove);
}

struct ant_t *initialize_data() {

  ant_t *local_ant = malloc_(sizeof(ant_t));
  local_ant->pheromones = malloc_(sizeof(double) * problem->nof_vertices);
  local_ant->phero_vars = malloc_(sizeof(double) * problem->nof_vertices);

  for (int i = 0; i < problem->nof_vertices; i++) {
    local_ant->pheromones[i] = malloc_(sizeof(double) * problem->nof_vertices);
    local_ant->phero_vars[i] = malloc_(sizeof(double) * problem->nof_vertices);
    for (int j = 0; j < problem->nof_vertices; j++) {
      local_ant->pheromones[i][j] = 0;
      local_ant->phero_vars[i][j] = 0;
      if (!problem->adj_matrix[i][j]) {
        local_ant->pheromones[i][j] = 1;
      }
    }
  }
  local_ant->best_ant = malloc_(sizeof(gcp_solution_t));
  local_ant->best_ant->color_of = malloc_(sizeof(int) * problem->nof_vertices);
  local_ant->best_ant->nof_confl_vertices = INT_MAX;
  local_ant->best_ant->nof_colors = problem->colors;
  local_ant->best_ant->spent_time_ls = 0;

  local_ant->best_colony = malloc_(sizeof(gcp_solution_t));
  local_ant->best_colony->color_of =
      malloc_(sizeof(int) * problem->nof_vertices);
  local_ant->best_colony->nof_confl_vertices = INT_MAX;
  local_ant->best_colony->nof_colors = problem->colors;
  local_ant->best_colony->spent_time_ls = 0;

  local_ant->ant_k = malloc_(sizeof(gcp_solution_t) * problem->nof_vertices);
  local_ant->ant_k->color_of = malloc_(sizeof(int) * problem->nof_vertices);
  local_ant->ant_k->nof_colors = problem->colors;
  local_ant->ant_k->spent_time_ls = 0;

  if (get_flag(problem->flags, FLAG_MEMORY)) {
    local_ant->ant_memory_remove = malloc_(sizeof(gcp_solution_t));
    local_ant->ant_memory_remove->color_of =
        malloc_(sizeof(int) * problem->nof_vertices);
    local_ant->ant_memory_remove->nof_colors = problem->colors;
    local_ant->ant_memory_remove->spent_time_ls = 0;

    local_ant->ant_memory_insert->color_of =
        malloc_(sizeof(int) * problem->nof_vertices);
    local_ant->ant_memory_insert->nof_colors = problem->colors;
    local_ant->ant_memory_insert->spent_time_ls = 0;
  }

  return local_ant;
}

void update_pheromone_trails_memory(gcp_solution_t *local_ant_memory_insert,
                                    gcp_solution_t *local_ant_memory_remove,
                                    double **pheromone) {
  int i, j;

  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {

      if (!problem->adj_matrix[i][j]) {
        if (local_ant_memory_insert->color_of[i] ==
            local_ant_memory_insert->color_of[j])
          pheromone[i][j] *= 1 + aco_info->delta;

        if ((get_flag(problem->flags, FLAG_MEMORY_REMOVE)) &&
            (local_ant_memory_remove->color_of[i] ==
             local_ant_memory_remove->color_of[j]))
          pheromone[i][j] *= 1 - aco_info->delta;
      }
    }
  }
}

void update_var_phero(gcp_solution_t *solution, double **phero_var) {

  int i, j;

  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {
      if (!problem->adj_matrix[i][j] &&
          (solution->color_of[i] == solution->color_of[j])) {
        phero_var[i][j] += (solution->nof_confl_vertices == 0)
                               ? 1
                               : 1.0 / solution->nof_confl_vertices;
      }
    }
  }
}

void update_pheromone_trails_scheme_1(gcp_solution_t *local_best_ant,
                                      gcp_solution_t *local_best_colony,
                                      double **pheromone, double **phero_var) {

  int i, j;

#if defined DEBUG
  fprintf(stderr, "update_pheromone_trails_scheme_1.\n");
#endif

  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {
      pheromone[i][j] *= aco_info->rho;
      pheromone[i][j] += phero_var[i][j];

      if (!problem->adj_matrix[i][j]) {
        if (local_best_ant->color_of[i] == local_best_ant->color_of[j])
          pheromone[i][j] += (local_best_ant->nof_confl_vertices == 0)
                                 ? 1
                                 : 1.0 / local_best_ant->nof_confl_vertices;
        if (local_best_colony->color_of[i] == local_best_colony->color_of[j])
          pheromone[i][j] += (local_best_colony->nof_confl_vertices == 0)
                                 ? 1
                                 : 1.0 / local_best_colony->nof_confl_vertices;
      }

      phero_var[i][j] = 0;
    }
  }
}

void update_pheromone_trails_scheme_2(gcp_solution_t *local_best_ant,
                                      gcp_solution_t *local_best_colony,
                                      double **pheromone) {
  for (int i = 0; i < problem->nof_vertices; i++) {
    for (int j = 0; j < problem->nof_vertices; j++) {
      pheromone[i][j] *= aco_info->rho;
      if (!problem->adj_matrix[i][j]) {
        if (local_best_ant->color_of[i] == local_best_ant->color_of[j]) {
          pheromone[i][j] += (local_best_ant->nof_confl_vertices == 0)
                                 ? 1
                                 : 1.0 / local_best_ant->nof_confl_vertices;
        }
        if (local_best_colony->color_of[i] == local_best_colony->color_of[j]) {
          pheromone[i][j] += (local_best_colony->nof_confl_vertices == 0)
                                 ? 1
                                 : 1.0 / local_best_colony->nof_confl_vertices;
        }
      }
    }
  }
}

void update_pheromone_trails_scheme_3(int cycle, gcp_solution_t *local_best_ant,
                                      gcp_solution_t *local_best_colony,
                                      double **pheromone) {

  int i, j;
  int gap = 0;

#if defined DEBUG
  fprintf(stderr, "update_pheromone_trails_scheme_3.\n");
#endif

  if (cycle % aco_info->gap == 0)
    gap = cycle / aco_info->gap;

  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = i; j < problem->nof_vertices; j++) {
      pheromone[i][j] *= aco_info->rho;
      if (!problem->adj_matrix[i][j]) {
        if (gap) {
          if (local_best_ant->color_of[i] == local_best_ant->color_of[j])
            pheromone[i][j] += (local_best_ant->nof_confl_vertices == 0)
                                   ? 1
                                   : 1.0 / local_best_ant->nof_confl_vertices;
        } else {
          if (local_best_colony->color_of[i] == local_best_colony->color_of[j])
            pheromone[i][j] +=
                (local_best_colony->nof_confl_vertices == 0)
                    ? 1
                    : 1.0 / local_best_colony->nof_confl_vertices;
        }
      }
      pheromone[j][i] = pheromone[i][j];
    }
  }
  gap--;
}

/* END Functions to help updating pheromone */

void construct_solutions(int cycle, ant_t **local_ant, ant_fixed_k_t **fixed_k,
                         tabucol_conflicts_t **tabucol_conflicts) {
  (*tabucol_conflicts) = malloc_(sizeof(tabucol_conflicts_t));
  (*tabucol_conflicts)->conf_position = NULL;
  (*tabucol_conflicts)->conflicts = NULL;
  (*tabucol_conflicts)->nodes_in_conflict = NULL;
  (*tabucol_conflicts)->tabu_status = NULL;
  aco_memory_t *memory = NULL;
  gcp_solution_t *ant_memory;
  int k;
  int ants = aco_info->ants;

  (*local_ant)->best_colony->f1 = INT_MAX * 1.0;
  for (k = 0; k < ants; k++) {

    (*local_ant)->ant_k->total_cycles = cycle;

    ant_fixed_k(fixed_k, (*local_ant)->ant_k, (*local_ant)->pheromones);

    /* Apply local search in all ants */
    if (get_flag(problem->flags, FLAG_TABUCOL_ALL_ANTS) &&
        ((*local_ant)->ant_k->nof_confl_vertices != 0) &&
        (tabucol_info->cycles > 0)) {
#if defined DEBUG
      // fprintf(stderr, "FLAG_TABUCOL_ALL_ANTS\n");
#endif

      tabucol(tabucol_conflicts, (*local_ant)->ant_k, tabucol_info->cycles,
              tabucol_info->tl_style);
    }

    (*local_ant)->ant_k->spent_time =
        current_time_secs(TIME_FINAL, time_initial);

    if ((*local_ant)->ant_k->f1 < (*local_ant)->best_colony->f1) {
      cpy_solution((*local_ant)->ant_k, (*local_ant)->best_colony);
      (*local_ant)->best_colony->cycles_to_best = cycle;
      (*local_ant)->best_colony->time_to_best = (*local_ant)->ant_k->spent_time;
    }

    if (aco_info->pheromone_scheme == PHEROMONE_SCHEME_1)
      update_var_phero((*local_ant)->ant_k, (*local_ant)->phero_vars);

    if (((*local_ant)->ant_k->nof_confl_vertices == 0) ||
        ((get_flag(problem->flags, FLAG_TIME)) &&
         (problem->time <= current_time_secs(TIME_FINAL, time_initial))))
      break;
  }

  if (!(get_flag(problem->flags, FLAG_TABUCOL_ALL_ANTS)) &&
      ((*local_ant)->best_colony->nof_confl_vertices != 0) &&
      (tabucol_info->cycles > 0)) {
#if defined DEBUG
    fprintf(stderr, "FLAG_TABUCOL_BEST_ANT\n");
#endif
    tabucol(tabucol_conflicts, (*local_ant)->best_colony, tabucol_info->cycles,
            tabucol_info->tl_style);
    (*local_ant)->best_colony->spent_time =
        current_time_secs(TIME_FINAL, time_initial);
  }

  if (get_flag(problem->flags, FLAG_MEMORY)) {
    ant_memory = malloc_(sizeof(gcp_solution_t));
    ant_memory->color_of = malloc_(sizeof(int) * problem->nof_vertices);
    cpy_solution((*local_ant)->best_colony, ant_memory);
    insert_into_memory(ant_memory, (*local_ant)->ant_memory_insert,
                       (*local_ant)->ant_memory_remove, memory);
  }
}

void execute_colorant(ant_t **local_ant, ant_fixed_k_t **ant_fixed_k,
                      tabucol_conflicts_t **tabucol_conflicts, int thread) {

  int cycle = 0;
  int converg = 0;
  int change = 0;
  int cycle_phero = 0;

  (*local_ant) = initialize_data();
  (*ant_fixed_k) = afk_initialize_data(aco_info->alpha, aco_info->beta);
  (*local_ant)->best_ant->stop_criterion = 0;

  while (!terminate_conditions((*local_ant)->best_ant, cycle, converg)) {

    cycle = cycle + 1;
    converg = converg + 1;
    cycle_phero = cycle_phero + 1;

    construct_solutions(cycle, local_ant, ant_fixed_k, tabucol_conflicts);

    if ((*local_ant)->best_colony->nof_confl_vertices <
        (*local_ant)->best_ant->nof_confl_vertices) {

      cpy_solution((*local_ant)->best_colony, (*local_ant)->best_ant);
      (*local_ant)->best_ant->cycles_to_best = cycle;
      (*local_ant)->best_ant->time_to_best =
          (*local_ant)->best_colony->spent_time;
      converg = 0;
      change = 1;
    }

    pthread_barrier_wait(&threads_barrier);

    // atualizar a melhor global
    pthread_mutex_lock(&global_best_ant_mutex);
    if (global_best_ant->nof_confl_vertices >
        (*local_ant)->best_ant->nof_confl_vertices) {
      cpy_solution((*local_ant)->best_ant, global_best_ant);
      global_best_ant->time_to_best = (*local_ant)->best_ant->spent_time;
    }
    pthread_mutex_unlock(&global_best_ant_mutex);

    pthread_barrier_wait(&threads_barrier);
    cpy_solution(global_best_ant, (*local_ant)->best_ant);

    switch (aco_info->pheromone_scheme) {
    case PHEROMONE_SCHEME_1:
      update_pheromone_trails_scheme_1(
          (*local_ant)->best_ant, (*local_ant)->best_colony,
          (*local_ant)->pheromones, (*local_ant)->phero_vars);
      break;
    case PHEROMONE_SCHEME_2:
      update_pheromone_trails_scheme_2((*local_ant)->best_ant,
                                       (*local_ant)->best_colony,
                                       (*local_ant)->pheromones);
      break;
    case PHEROMONE_SCHEME_3:
      update_pheromone_trails_scheme_3(cycle, (*local_ant)->best_ant,
                                       (*local_ant)->best_colony,
                                       (*local_ant)->pheromones);
      break;
    }

    if (get_flag(problem->flags, FLAG_MEMORY)) {
      update_pheromone_trails_memory((*local_ant)->ant_memory_insert,
                                     (*local_ant)->ant_memory_remove,
                                     (*local_ant)->pheromones);
    }

    if (get_flag(problem->flags, FLAG_VERBOSE)) {
      fprintf(problem->fileout,
              "Cycle %d - Conflicts found: %d (edges), %d (vertices)\n", cycle,
              (*local_ant)->best_ant->nof_confl_edges,
              (*local_ant)->best_ant->nof_confl_vertices);
      fflush(stdout);
    }

    if ((*local_ant)->best_ant->nof_confl_vertices == 0) {
      (*local_ant)->best_ant->spent_time =
          current_time_secs(TIME_FINAL, time_initial);
      (*local_ant)->best_ant->total_cycles = cycle;
      (*local_ant)->best_ant->stop_criterion = STOP_BEST;
      break;
    }

    if (get_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA) &&
        ((cycle % aco_info->iterations_alpha_beta) == 0)) {

      aco_info->gamma = change ? (1 - aco_info->omega) * aco_info->gamma
                               : (1 + aco_info->omega) * aco_info->gamma;

      aco_info->alpha = aco_info->alpha_base * aco_info->gamma;
      aco_info->beta = aco_info->beta_base * (1 - aco_info->gamma);

      if ((aco_info->alpha <= 0) && (aco_info->beta > 0))
        aco_info->gamma = (1 + aco_info->omega) * aco_info->gamma;
      else if ((aco_info->alpha > 0) && (aco_info->beta <= 0))
        aco_info->gamma = (1 - aco_info->omega) * aco_info->gamma;

      aco_info->alpha = aco_info->alpha_base * aco_info->gamma;
      aco_info->beta = aco_info->beta_base * (1 - aco_info->gamma);

      change = 0;
    }

    if ((get_flag(problem->flags, FLAG_DIFF_TABUCOL_SCHEME)) &&
        ((cycle % tabucol_info->diff_scheme_iterations) == 0)) {
      tabucol_info->tl_style = tabucol_info->tl_style == TABUCOL_REACTIVE
                                   ? TABUCOL_DYNAMIC
                                   : TABUCOL_REACTIVE;
    }
    if ((get_flag(problem->flags, FLAG_CHANGE_PHEROMONE_SCHEME)) &&
        (converg >= aco_info->change_phero_scheme_iterations) &&
        (cycle_phero >= aco_info->change_phero_scheme_iterations)) {

      aco_info->pheromone_scheme++;
      if (aco_info->pheromone_scheme > 3)
        aco_info->pheromone_scheme = 1;

      cycle_phero = 0;
    }
  }

  return;
}
