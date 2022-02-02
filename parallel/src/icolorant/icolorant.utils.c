
/*************************************************************************
 * Created: Qua 27 Jan 2011 15:17:39 BRST
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 * ModifiedAt: 02/02/2022 - 06:54
 * Modified By: Thiao Issao Yasunaka, thiagoyasunaka@hotmail.com
 *
 *************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "icolorant.h"
#include "../tabucol.h"
#include "../ant_fixed_k.h"
#include "../util.h"

void colorant_printbanner(void) {

  char *schemes[] = {"All ants + Best ant + Best colony",
                     "Best ant + Best colony", "Best ant + Best colony (gap)"};

  fprintf(problem->fileout, "COLORANT\n");
  fprintf(problem->fileout,
          "-------------------------------------------------\n");
  fprintf(problem->fileout, "Parameters:\n");
  fprintf(problem->fileout, "  Alpha..................................: %.2f\n",
          aco_info->alpha);
  fprintf(problem->fileout, "  Beta...................................: %.2f\n",
          aco_info->beta);
  fprintf(problem->fileout, "  Rho....................................: %.2f\n",
          aco_info->rho);
  if (!(get_flag(problem->flags, FLAG_ANTS_RATIO)))
    fprintf(problem->fileout, "  Ants...................................: %i\n",
            aco_info->ants);
  else
    fprintf(
        problem->fileout,
        "  Ants...................................: %i (%i of %i - vertices)\n",
        aco_info->ants, aco_info->ratio, problem->nof_vertices);
  fprintf(problem->fileout, "  Pheromone scheme.......................: %s\n",
          schemes[aco_info->pheromone_scheme - 1]);

  if (get_flag(problem->flags, FLAG_CHANGE_PHEROMONE_SCHEME))
    fprintf(problem->fileout,
            "  Change pheromone scheme after %i iterations.\n",
            aco_info->change_phero_scheme_iterations);

  if (get_flag(problem->flags, FLAG_MEMORY)) {
    if (!(get_flag(problem->flags, FLAG_MEMORY_RATIO))) {
      fprintf(problem->fileout,
              "  Memory Usage:\n\tMemory size......................: %i\n",
              aco_info->memory_size);
    } else {
      fprintf(problem->fileout,
              "  Memory Usage:\n\tMemory size......................: %i (%i of "
              "%i - ants)\n",
              aco_info->memory_size, aco_info->memory_ratio, aco_info->ants);
    }
    fprintf(problem->fileout, "\tDelta............................: %.2f\n",
            aco_info->delta);
  }

  if (aco_info->pheromone_scheme == PHEROMONE_SCHEME_3)
    fprintf(problem->fileout,
            "  Pheromone Scheme 3:\n\tGap..............................: %i\n",
            aco_info->gap);

  if (get_flag(problem->flags, FLAG_REUSE_COLOR)) {
    fprintf(
        problem->fileout,
        "  Try to reuse colors:\n\tX................................: %.2f\n",
        aco_info->x);
    fprintf(problem->fileout, "\tY................................: %.2f\n",
            aco_info->y);
  }

  if (get_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA)) {
    fprintf(
        problem->fileout,
        "  Change alpha and beta:\n\tGamma............................: %.2f\n",
        aco_info->gamma);
    fprintf(problem->fileout, "\tOmega............................: %.2f\n",
            aco_info->omega);
    fprintf(problem->fileout, "\tChange alpha and beta after %i iterations.\n",
            aco_info->iterations_alpha_beta);
  }
}

void colorant_show_solution(void) {
  if (get_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA)) {
    fprintf(problem->fileout,
            "-------------------------------------------------\n");
    fprintf(problem->fileout, "Alpha..............................: %.2f\n",
            aco_info->alpha);
    fprintf(problem->fileout, "Beta...............................: %.2f\n",
            aco_info->beta);
    fprintf(problem->fileout, "Rho................................: %.2f\n",
            aco_info->rho);
  }
}

void print_memory(aco_memory_t *memory, gcp_solution_t *ant_memory_remove) {
  aco_memory_t *lmemory = memory;
  int v, count = 1;
  for (; lmemory; lmemory = lmemory->tail) {
    fprintf(problem->fileout, "Item: %i\n", count++);
    fprintf(problem->fileout, "No. of conflicting edges: %d\n",
            lmemory->head->nof_confl_edges);
    fprintf(problem->fileout, "No. of    conflicting vertices: %d\n",
            lmemory->head->nof_confl_vertices);
    fprintf(problem->fileout, "Color:\n");
    for (v = 0; v < problem->nof_vertices; v++)
      fprintf(problem->fileout, "%i, ", lmemory->head->color_of[v]);
    fprintf(problem->fileout, "\n\n");
  }

  fprintf(problem->fileout, "Removed:\n");
  fprintf(problem->fileout, "No. of conflicting edges: %d\n",
          ant_memory_remove->nof_confl_edges);
  fprintf(problem->fileout, "No. of  conflicting vertices: %d\n",
          ant_memory_remove->nof_confl_vertices);
  fprintf(problem->fileout, "Color:\n");
  for (v = 0; v < problem->nof_vertices; v++)
    fprintf(problem->fileout, "%i, ", ant_memory_remove->color_of[v]);
  fprintf(problem->fileout, "\n\n");
}