/***********************************************************
 * Created: Seg 10 Out 2011 17:17:12 BRT
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 * ANT_FIXED_K
 * * Constructive method for each ant in an ACO algorithm for k-GCP
 *
 ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "icolor.h"
#include "./icolorant/icolorant.h"
#include "ant_fixed_k.h"
#include "util.h"

/* ANT_FIXED_K data */
static double *probb;
static double **trail;
static float alpha;
static float beta;

static int choose_vertex(
    int neighbors_by_color[problem->nof_vertices][problem->colors + 1],
    int *color_of) {

  int v = 0, i, dsat, maxdsat;

  /* Choose vertex with maximum degree of saturation */
  maxdsat = -1;
  for (i = 0; i < problem->nof_vertices; i++) {
    if (color_of[i] == -1) {
      dsat = neighbors_by_color[i][problem->colors];
      if (dsat > maxdsat) {
        v = i;
        maxdsat = dsat;
      }
    }
  }
  return v;
}

static void calculate_probbs(
    int v, int *color_of, int size_color[problem->colors],
    int neighbors_by_color[problem->nof_vertices][problem->colors + 1]) {

  int c;
  double sum, traill, totalsum, neighbors;

  totalsum = 0;

  for (c = 0; c < problem->colors; c++) {

    probb[c] = 0;

    sum = trail[c][v];

    if (get_flag(problem->flags, FLAG_REUSE_COLOR)) {
#if defined DEBUG
      fprintf(stderr, "trying to reuse color.\n");
#endif
      if (size_color[c] == 0) {
        traill = aco_info->y;
      } else {
        if (neighbors_by_color[v][c] == 0) {
          traill = aco_info->x;
        } else {
          traill = sum / size_color[c];
        }
      }
    } else
      traill = (size_color[c] == 0) ? 1 : sum / size_color[c];

    neighbors = neighbors_by_color[v][c] + 1;
    neighbors = 1.0 / neighbors;

    probb[c] = pow(traill, alpha) * pow(neighbors, beta);

    totalsum += probb[c];
  }

  totalsum = (totalsum == 0) ? 1 : totalsum;

  /* To avoid a new 'for'  */
  probb[problem->colors] = totalsum;
}

static int choose_color(void) {

  int i;
  double p, last, div;

  div = probb[problem->colors];

#if defined LRAND
  // p = (double) RANDOM_UNIT() / INT_MAX;
  RANDOM_UNIT(problem->buffer, p, double);
  p = p / INT_MAX;
#elif defined NRAND
  // p = (double) RANDOM_UNIT(problem->seed) / INT_MAX;
  RANDOM_UNIT(problem->seed, problem->buffer, p, double);
  p = p / INT_MAX;
#endif

  last = 0;
  for (i = 0; i < problem->colors; i++) {
    last += (probb[i] / div);
    if (p <= last) {
      return i;
    }
  }
  /* When it reaches here, it means that p == 1 */
  return problem->colors - 1;
}

void afk_initialize_data(float p_alpha, float p_beta) {

  int i, j;

  probb = malloc_(sizeof(double) * (problem->colors + 1));
  trail = malloc_(sizeof(double *) * problem->colors);

  for (i = 0; i < problem->colors; i++) {
    trail[i] = malloc_(sizeof(double) * problem->nof_vertices);
    for (j = 0; j < problem->nof_vertices; j++) {
      trail[i][j] = 0;
    }
  }

  alpha = p_alpha;
  beta = p_beta;
}

void ant_fixed_k(gcp_solution_t *solution, double **pheromone) {

  int i, j;
  int color = 0;   /* number of colors to be used */
  int colored = 0; /* number of colored vertex */
  int v;           /* vertex to be colored */

  int confl_vertices[problem->nof_vertices];
  int neighbors_by_color[problem->nof_vertices][problem->colors + 1];
  int size_color[problem->colors];

  double sum = 0.0;

  solution->nof_colors = problem->colors;

  /* Initializing auxiliary arrays */
  for (i = 0; i < problem->nof_vertices; i++) {
    solution->color_of[i] = -1;
    size_color[i] = 0;
    confl_vertices[i] = 0;
    for (j = 0; j < problem->colors; j++) {
      neighbors_by_color[i][j] = 0;
      trail[j][i] = 0;
    }
    neighbors_by_color[i][problem->colors] = 0;
  }

  solution->nof_confl_edges = 0;
  solution->nof_confl_vertices = 0;

  while (colored < problem->nof_vertices) {

    /* Chose a vertex to be colored */
    v = choose_vertex(neighbors_by_color, solution->color_of);

    /* Calculate colors probabilities */
    calculate_probbs(v, solution->color_of, size_color, neighbors_by_color);

    /* Choose a color to be assigned to v */
    color = choose_color();

    solution->color_of[v] = color;
    size_color[color]++;
    colored++;

    /* Update informations about conflicts and saturation degree */
    int conf = solution->nof_confl_edges;
    for (i = 0; i < problem->nof_vertices; i++) {

      /* trail keeps the pheromone between a vertex and all the vertex
       * already colored with each color */
      trail[color][i] += pheromone[v][i];

      if (problem->adj_matrix[v][i]) {
        /* update degree of saturation: */
        if (neighbors_by_color[i][color] == 0) {
          neighbors_by_color[i][problem->colors]++;
        }
        /* now <i> has a neighbor colored with <color> */
        neighbors_by_color[i][color]++;

        /* if a neighbor of <v> is colored with <color>, there is a
         * conflicting edge between them */
        if (solution->color_of[i] == color) {
          solution->nof_confl_edges++;
          if (confl_vertices[i] == 0) {
            confl_vertices[i] = 1;
            solution->nof_confl_vertices++;
          }
        }
      }
    }
    /* if any new conflicting edge was created, <v> is a conflicting
     * vertex */
    if (conf != solution->nof_confl_edges) {
      if (confl_vertices[v] == 0) {
        confl_vertices[v] = 1;
        solution->nof_confl_vertices++;
      }
    }
  }

  for (i = 0; i < problem->nof_vertices; i++)
    sum += confl_vertices[i] * (1.0 / problem->degree[i]);

  solution->h1 = sum * (1.0 / problem->nof_vertices);
  solution->f1 = solution->nof_confl_vertices - solution->h1;

  solution->spent_time = current_time_secs(TIME_FINAL, time_initial);
}
